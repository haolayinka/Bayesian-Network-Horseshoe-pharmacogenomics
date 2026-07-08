#!/usr/bin/env Rscript
#
# gibbs_sampler.R
#
# Implements the Gibbs sampler specified in Section 5 of pathway.tex for
# the network-structured horseshoe drug-response model (Section 4):
#
#   y_id ~ Normal(mu_id, sigma_d^2),   mu_id = alpha_{c_i,d} + sum_g m_ig beta_gd
#   beta_gd ~ Normal(beta_tilde_gd, eta_gd^2)                [optional group layer]
#   beta_tilde_.d ~ Normal(0, kappa_d^2 (D_W - delta W)^{-1})  [GMRF network prior]
#   beta_tilde_gd ~ Normal(0, lambda_gd^2 zeta_d^2)            [horseshoe]
#   lambda_gd ~ Half-Cauchy(0, 1)                              [LOCAL, via IG augmentation]
#   1/kappa_d^2 ~ Gamma(a, b),  1/zeta_d^2 ~ Gamma(a, b)        [GLOBAL/GROUP, GIGG-style]
#
# Every step below is a direct conjugate draw (Normal, Inverse-Gamma, or
# Gamma); there is no Metropolis-Hastings step anywhere in this sampler.
# This file defines the sampler as a set of functions, intended to be
# sourced by a driver script (run_simulation_test.R or run_real_data.R)
# rather than run directly.

suppressWarnings(suppressMessages({
  library(Matrix)
}))

# ----------------------------------------------------------------------- #
# Inverse-Gamma sampling helper (R has no built-in rinvgamma without
# extra packages; implemented directly from rgamma to avoid an additional
# dependency, since IG(a, b) <=> 1 / Gamma(a, rate = b)). Still used for
# the LOCAL horseshoe shrinkage parameter lambda_gd, which keeps its
# original Inverse-Gamma structure (see header note above).
# ----------------------------------------------------------------------- #

rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

# ----------------------------------------------------------------------- #
# gibbs_sampler()
#
# Arguments:
#   y          : N x D matrix of response values (y[i, d]); NA allowed for
#                missing (model, drug) pairs (e.g. if not every drug was
#                screened on every cell line) -- NA entries are skipped in
#                every likelihood-related sum via the obs_mask below.
#   M          : N x G binary mutation matrix (m_ig), genes as columns.
#   group      : length-N integer vector of group membership c_i in
#                1..K (e.g. tissue type); set to rep(1, N) to disable the
#                optional group-pooling layer (every sample falls in one
#                group, equivalent to no pooling).
#   W          : G' x G' symmetric, zero-diagonal pathway adjacency matrix.
#                Must use a subset of the gene names in colnames(M); genes
#                in M but not in W are still modeled, just without a
#                network-correlated prior (see beta_tilde_idx below).
#   delta      : fixed GMRF smoothing parameter (scalar in (0,1)).
#   use_group_layer, use_network_layer, use_horseshoe_layer : booleans
#                controlling which of the three optional layers are
#                active. All default to TRUE except use_group_layer
#                (default TRUE for backward compatibility, but FALSE --
#                i.e. no group pooling -- is the data-recommended setting
#                for the real dataset; see pathway.tex Section 4 and the
#                project progress log). Setting use_network_layer = FALSE
#                or use_horseshoe_layer = FALSE fits the corresponding
#                ablation comparator specified in pathway.tex Phase 3.
#   flat_prior_var : prior variance for beta_tilde when
#                use_horseshoe_layer = FALSE (default 100).
#   n_iter     : total number of Gibbs iterations.
#   burn_in    : number of initial iterations discarded before returning
#                samples.
#   thin       : keep every `thin`-th post-burn-in iteration.
#   verbose    : print progress every `verbose` iterations (0 = silent).
#   seed       : RNG seed for reproducibility.
#   gigg_a, gigg_b : Gamma-prior hyperparameters for the precision of
#                kappa2 and zeta2 (GIGG-style fix; see the file header
#                MODEL CORRECTION note).
#
# Returns a list of posterior sample arrays (post burn-in, thinned):
#   alpha       : [K, D, n_kept]
#   beta_tilde  : [G, D, n_kept]
#   sigma2      : [D, n_kept]
#   kappa2      : [D, n_kept]  (all-NA for a kept iteration if use_network_layer=FALSE)
#   lambda2     : [G, D, n_kept]  (all-NA if use_horseshoe_layer=FALSE)
#   zeta2       : [D, n_kept]  (all-NA if use_horseshoe_layer=FALSE)
#   gene_names, group_levels, keep_iters
#   config      : list recording exactly which layers/hyperparameters
#                produced this fit, so a saved fit is self-documenting
#                when comparing multiple ablation configurations later.
# ----------------------------------------------------------------------- #

gibbs_sampler <- function(y, M, group, W, delta = 0.8,
                           use_group_layer = TRUE,
                           use_network_layer = TRUE,
                           use_horseshoe_layer = TRUE,
                           flat_prior_var = 100,
                           n_iter = 4000, burn_in = 2000, thin = 2,
                           verbose = 200, seed = 1,
                           gigg_a = 1, gigg_b = 1) {
  # gigg_a, gigg_b: shape/rate hyperparameters for the Gamma prior on the
  # precision of each GLOBAL shrinkage parameter (1/kappa2_d, 1/zeta2_d).
  # Gamma(1,1) is the default used during validation on real data (see
  # project progress log); both parameters stabilized at sensible,
  # non-degenerate values under this choice across all genes and drugs
  # tested. As with delta, a small sensitivity grid over (gigg_a, gigg_b)
  # is recommended before treating any single choice as final for the
  # paper's results.
  #
  # use_network_layer, use_horseshoe_layer: set either to FALSE to fit the
  # corresponding ABLATION comparator specified in pathway.tex Phase 3
  # (Section 3.4): a no-network variant (genes exchangeable, closer to
  # Samorodnitsky et al. 2020's pooling structure) or a no-horseshoe
  # variant (flat priors on gene effects, no sparsity layer). Both
  # default to TRUE (the full model); setting both FALSE simultaneously
  # is allowed and fits beta_tilde with ONLY the flat prior below -- a
  # useful third comparator isolating the data term alone. When a layer
  # is disabled, its associated parameters (kappa2/W/delta for the
  # network layer; lambda2/zeta2/nu for the horseshoe layer) are not
  # estimated and are excluded from the returned output for that fit,
  # rather than being computed and silently ignored.
  #
  # flat_prior_var: the prior variance used for beta_tilde when
  # use_horseshoe_layer = FALSE, i.e. beta_tilde_gd ~ Normal(0,
  # flat_prior_var) independently per gene. Default 100, matching the
  # flat prior already used for alpha_kd elsewhere in this model -- a
  # genuinely uninformative prior across mutation-effect coefficients'
  # realistic range (which is empirically O(0.01) to O(1), per the real
  # gene/coefficient values observed during validation), not a
  # tightly-tuned value chosen to make this comparator look artificially
  # good or bad.

  set.seed(seed)

  N <- nrow(y)
  D <- ncol(y)
  G <- ncol(M)
  gene_names <- colnames(M)

  stopifnot(nrow(M) == N)
  stopifnot(length(group) == N)

  K <- length(unique(group))
  group_idx <- as.integer(factor(group))  # remap to 1..K contiguous

  obs_mask <- !is.na(y)  # N x D logical: TRUE where y is observed

  # --- Network prior setup ---
  # Genes in M that also appear in W get the GMRF network prior on
  # beta_tilde; genes in M absent from W (e.g. the 83 unmapped genes in
  # the real panel) get horseshoe-only shrinkage, handled by giving them
  # an "infinite-variance" (i.e. zero-precision contribution) row/column
  # in the network precision matrix below rather than excluding them from
  # the model entirely. If use_network_layer = FALSE (the no-network
  # ablation comparator, pathway.tex Phase 3), G_net is forced to 0
  # regardless of W, which disables the network term entirely through the
  # same `if (G_net > 0)` gate already used for genes absent from W --
  # no separate code path is needed for this ablation.
  if (use_network_layer) {
    net_genes <- intersect(gene_names, rownames(W))
  } else {
    net_genes <- character(0)
  }
  net_idx <- match(net_genes, gene_names)         # positions within 1..G
  G_net <- length(net_genes)

  if (G_net > 0) {
    W_sub <- W[net_genes, net_genes, drop = FALSE]
    D_W <- diag(rowSums(W_sub))
    Prec_net_sub <- D_W - delta * W_sub  # (D_W - delta*W), to be scaled by 1/kappa_d^2
    # Verify positive-definiteness once, up front, rather than discovering
    # a failure deep inside the sampler loop -- this mirrors the check
    # already performed during network construction (see the project
    # progress log) and is repeated here defensively in case W or delta
    # changes between sessions.
    eigvals_check <- eigen(Prec_net_sub, symmetric = TRUE, only.values = TRUE)$values
    if (min(eigvals_check) <= 1e-8) {
      stop(sprintf(
        "D_W - delta*W is not positive definite at delta=%.3f (min eigenvalue = %.6e). Re-run the network construction step or choose a different delta from the verified grid {0.5, 0.7, 0.8, 0.9}.",
        delta, min(eigvals_check)
      ))
    }
  } else {
    Prec_net_sub <- matrix(numeric(0), 0, 0)
  }

  # --- Initialize parameters ---
  alpha <- matrix(0, K, D)
  beta_tilde <- matrix(0, G, D)
  beta <- array(0, dim = c(K, G, D))  # group-level beta (used only if use_group_layer)
  sigma2 <- rep(1, D)
  eta2 <- matrix(1, G, D)
  kappa2 <- rep(1, D)
  lambda2 <- matrix(1, G, D)
  zeta2 <- rep(1, D)
  nu <- matrix(1, G, D)    # horseshoe LOCAL auxiliary variable (lambda2 stays Inverse-Gamma)
  # Note: the horseshoe GLOBAL auxiliary variable "xi" used in the original
  # Makalic-Schmidt half-Cauchy-via-IG-mixture scheme is no longer needed --
  # zeta2's precision is now drawn directly from a Gamma distribution (see
  # Step 7 below and the file header), which is already a complete
  # conjugate family requiring no further auxiliary-variable augmentation.

  # --- Storage for kept samples ---
  keep_iters <- seq(burn_in + 1, n_iter, by = thin)
  n_kept <- length(keep_iters)
  out_alpha <- array(NA, dim = c(K, D, n_kept))
  out_beta_tilde <- array(NA, dim = c(G, D, n_kept))
  out_sigma2 <- matrix(NA, D, n_kept)
  out_kappa2 <- matrix(NA, D, n_kept)
  out_lambda2 <- array(NA, dim = c(G, D, n_kept))
  out_zeta2 <- matrix(NA, D, n_kept)
  keep_pos <- 0L

  # Precompute, per drug, which rows (samples) are observed -- avoids
  # repeating the is.na() scan every iteration.
  obs_rows_by_drug <- lapply(seq_len(D), function(d) which(obs_mask[, d]))

  # ----------------------------------------------------------------- #
  # PERFORMANCE: precompute every per-drug quantity that does NOT depend
  # on the current iteration's sampled parameters, once, before the
  # iteration loop -- rather than recomputing it fresh on every one of
  # n_iter iterations. The dominant cost identified by profiling on real
  # data (see the project progress log) was crossprod(M_d) in the
  # no-group-layer path, which is a function of the fixed mutation data
  # alone and is therefore identical on every iteration for a given drug;
  # recomputing a ~951x219 crossprod thousands of times per drug was
  # responsible for roughly 60% of total runtime. M_d itself (the row
  # slice of M for drug d's observed samples) is precomputed for the same
  # reason: slicing is cheap per call but still unnecessary repeated work
  # across thousands of iterations.
  # ----------------------------------------------------------------- #
  M_d_list <- vector("list", D)
  y_d_list <- vector("list", D)
  grp_d_list <- vector("list", D)
  data_prec_full_unscaled_list <- vector("list", D)  # crossprod(M_d), no-group-layer path only

  for (d in seq_len(D)) {
    rows_d <- obs_rows_by_drug[[d]]
    if (length(rows_d) == 0) next
    M_d_list[[d]] <- M[rows_d, , drop = FALSE]
    y_d_list[[d]] <- y[rows_d, d]
    grp_d_list[[d]] <- group_idx[rows_d]
    if (!use_group_layer) {
      data_prec_full_unscaled_list[[d]] <- crossprod(M_d_list[[d]])
    }
  }

  for (t in seq_len(n_iter)) {

    for (d in seq_len(D)) {
      rows_d <- obs_rows_by_drug[[d]]
      if (length(rows_d) == 0) next  # drug screened on zero retained samples; skip

      y_d <- y_d_list[[d]]
      M_d <- M_d_list[[d]]
      grp_d <- grp_d_list[[d]]

      # current fitted mutation contribution, used to back out the
      # intercept's residual and vice versa
      mu_mut_d <- M_d %*% beta_tilde[, d]  # using beta_tilde directly if no group layer;
                                            # overwritten below if group layer active
      if (use_group_layer) {
        beta_d_per_row <- beta[grp_d, , d, drop = FALSE]
        # recompute mu_mut_d using the group-specific beta for each row
        mu_mut_d <- rowSums(M_d * beta[grp_d, , d])
      }

      # ---------------------------------------------------------------- #
      # Step 1: intercepts alpha_{kd}
      # ---------------------------------------------------------------- #
      for (k in seq_len(K)) {
        rows_k <- rows_d[grp_d == k]
        if (length(rows_k) == 0) next
        idx_in_d <- which(grp_d == k)
        resid_k <- y_d[idx_in_d] - mu_mut_d[idx_in_d]
        prec_k <- length(rows_k) / sigma2[d] + 1 / 100
        mean_k <- (sum(resid_k) / sigma2[d]) / prec_k
        alpha[k, d] <- rnorm(1, mean_k, sqrt(1 / prec_k))
      }

      # recompute residual against updated intercepts
      resid_d <- y_d - alpha[grp_d, d]

      # ---------------------------------------------------------------- #
      # Step 2 & 3: group-level beta and group heterogeneity variance
      # (optional layer only)
      # ---------------------------------------------------------------- #
      if (use_group_layer) {
        for (k in seq_len(K)) {
          idx_in_d <- which(grp_d == k)
          if (length(idx_in_d) == 0) {
            beta[k, , d] <- rnorm(G, beta_tilde[, d], sqrt(eta2[, d]))
            next
          }
          M_k <- M_d[idx_in_d, , drop = FALSE]
          y_k <- resid_d[idx_in_d]
          # Conjugate Normal update for a G-vector regression coefficient
          # with independent Normal(beta_tilde, eta2) prior per gene.
          prior_prec <- diag(1 / eta2[, d], G, G)
          A_k <- crossprod(M_k) / sigma2[d] + prior_prec
          b_k <- crossprod(M_k, y_k) / sigma2[d] + (beta_tilde[, d] / eta2[, d])
          # Single Cholesky factorization serves both the posterior mean
          # and the random draw -- see the identical optimization and its
          # rationale documented above Step 4's B_d update.
          R_k <- chol(A_k)
          beta_mean <- backsolve(R_k, forwardsolve(t(R_k), b_k))
          beta[k, , d] <- as.numeric(
            beta_mean + backsolve(R_k, rnorm(G))
          )
        }
        # Step 3: group heterogeneity variance eta2_gd
        for (g in seq_len(G)) {
          ss <- sum((beta[, g, d] - beta_tilde[g, d])^2)
          eta2[g, d] <- rinvgamma(1, shape = K / 2 + 0.01, rate = 0.01 + ss / 2)
        }
      }

      # ---------------------------------------------------------------- #
      # Step 4: network-and-horseshoe-structured mean effects beta_tilde_.d
      # ---------------------------------------------------------------- #
      if (use_group_layer) {
        # Data-driven term comes from pooling across groups (sum over k of
        # beta_{kgd}), precision diag(1/eta2_gd) * K in the simplest
        # complete-pooling case; here we use the exact per-group precision
        # sum so unequal eta2_gd values, if they arise, are respected.
        data_prec_diag <- rep(0, G)
        data_b <- rep(0, G)
        for (k in seq_len(K)) {
          data_prec_diag <- data_prec_diag + 1 / eta2[, d]
          data_b <- data_b + beta[k, , d] / eta2[, d]
        }
      } else {
        # No-group variant: beta_tilde regresses directly on the data.
        # data_prec_full_unscaled_list[[d]] = crossprod(M_d), precomputed
        # once before the iteration loop since it does not depend on any
        # sampled parameter -- this was the single largest cost in
        # profiling (see comment above the precomputation block).
        data_prec_full <- data_prec_full_unscaled_list[[d]] / sigma2[d]
        data_b_full <- crossprod(M_d, resid_d) / sigma2[d]
      }

      # Horseshoe precision (diagonal): 1 / (lambda_gd^2 * zeta_d^2).
      # If use_horseshoe_layer = FALSE (the no-horseshoe ablation
      # comparator, pathway.tex Phase 3), beta_tilde instead gets a flat,
      # fixed Normal(0, flat_prior_var) prior per gene -- a genuinely
      # uninformative prior (no sparsity, no shrinkage), not a
      # disguised replacement for the horseshoe. lambda2/zeta2/nu are not
      # estimated in this case (Steps 6-7 are skipped below).
      if (use_horseshoe_layer) {
        hs_prec_diag <- 1 / (lambda2[, d] * zeta2[d])
      } else {
        hs_prec_diag <- rep(1 / flat_prior_var, G)
      }

      # Network precision (only on the G_net genes that map to a pathway):
      # full G x G precision matrix, zero outside the net_idx submatrix,
      # so genes with no pathway membership get zero contribution here
      # (they rely entirely on the horseshoe term above, as documented in
      # pathway.tex Section 4). If use_network_layer = FALSE (the
      # no-network ablation comparator), G_net is forced to 0 above, so
      # this matrix is always all-zero here and contributes nothing --
      # kappa2 is correspondingly not estimated (Step 5 is skipped below).
      net_prec_full <- matrix(0, G, G)
      if (G_net > 0) {
        net_prec_full[net_idx, net_idx] <- Prec_net_sub / kappa2[d]
      }

      if (use_group_layer) {
        B_d <- diag(data_prec_diag, G, G) + net_prec_full + diag(hs_prec_diag, G, G)
        b_d <- data_b
      } else {
        B_d <- data_prec_full + net_prec_full + diag(hs_prec_diag, G, G)
        b_d <- as.numeric(data_b_full)
      }

      # PERFORMANCE: B_d is symmetric positive-definite (it's a precision
      # matrix: sum of a data-precision term, the GMRF network precision,
      # and the horseshoe precision, each PSD or PD), so it admits a
      # single Cholesky factorization R (upper-triangular, B_d = R'R) that
      # serves BOTH purposes below -- solving for the posterior mean via
      # two triangular solves, and drawing the random Normal perturbation
      # via a single triangular solve against a standard Normal vector.
      # This replaces the previous solve(B_d) (a full LU-based matrix
      # inverse) followed by a SECOND, separate chol() of that inverse,
      # which was redundant: both decompositions describe the same
      # underlying covariance structure, so doing it once is mathematically
      # equivalent and meaningfully faster (this and the analogous Step 2
      # update were the dominant remaining cost after the crossprod fix
      # documented above, per profiling on real data).
      R_d <- chol(B_d)  # upper-triangular, B_d = t(R_d) %*% R_d
      beta_tilde_mean <- backsolve(R_d, forwardsolve(t(R_d), b_d))
      beta_tilde[, d] <- as.numeric(
        beta_tilde_mean + backsolve(R_d, rnorm(G))
      )

      # ---------------------------------------------------------------- #
      # Step 5: network smoothness variance kappa2_d
      #
      # GIGG-style fix (see file header): kappa2_d's precision
      # (1/kappa2_d) is drawn from a Gamma distribution rather than
      # kappa2_d itself from an Inverse-Gamma -- this breaks the
      # distributional symmetry with zeta2_d (Step 7) that was previously
      # causing both to collapse toward zero together. Derivation: with
      # beta_tilde_net | tau_kappa ~ N(0, Prec_net^{-1} / tau_kappa), the
      # likelihood contribution for tau_kappa = 1/kappa2_d is
      # Gamma(shape=G_net/2, rate=quad_form/2); combined with a
      # Gamma(gigg_a, gigg_b) prior, the conjugate posterior is
      # Gamma(gigg_a + G_net/2, gigg_b + quad_form/2). Verified numerically
      # against a directly grid-evaluated posterior before adoption.
      # ---------------------------------------------------------------- #
      if (G_net > 0) {
        bt_net <- beta_tilde[net_idx, d]
        quad_form <- as.numeric(t(bt_net) %*% Prec_net_sub %*% bt_net)
        tau_kappa <- rgamma(1, shape = gigg_a + G_net / 2, rate = gigg_b + quad_form / 2)
        kappa2[d] <- 1 / tau_kappa
      }

      # ---------------------------------------------------------------- #
      # Step 6: horseshoe local shrinkage (per gene)
      #
      # Unchanged: lambda_gd keeps its original Inverse-Gamma (half-Cauchy
      # via Makalic-Schmidt augmentation) structure. GIGG prescribes
      # leaving the LOCAL shrinkage parameter Inverse-Gamma (heavy-tailed,
      # lets individual genes escape shrinkage when supported by data) and
      # only changing the GLOBAL/GROUP-level parameter's family -- see
      # Step 7 below and the file header. Skipped entirely if
      # use_horseshoe_layer = FALSE (the no-horseshoe ablation comparator,
      # pathway.tex Phase 3): lambda2/nu are not estimated in that case,
      # since beta_tilde instead uses the flat prior set up in Step 4.
      # ---------------------------------------------------------------- #
      if (use_horseshoe_layer) {
        for (g in seq_len(G)) {
          lambda2[g, d] <- rinvgamma(1, shape = 1,
                                      rate = 1 / nu[g, d] + beta_tilde[g, d]^2 / (2 * zeta2[d]))
          nu[g, d] <- rinvgamma(1, shape = 1, rate = 1 + 1 / lambda2[g, d])
        }
      }

      # ---------------------------------------------------------------- #
      # Step 7: horseshoe global shrinkage (per drug)
      #
      # GIGG-style fix (see file header and Step 5 above, applied
      # identically here): zeta2_d's precision is drawn from a Gamma
      # distribution instead of zeta2_d itself from an Inverse-Gamma.
      # zeta2_d has shape=(G+1)/2 (G=219 in the real panel, i.e. shape~110)
      # -- an even more severe instance of the Step 5 pathology, confirmed
      # on real data: zeta2_d collapsed to exactly 0 even after kappa2_d
      # was stabilized by the Step 5 fix alone, until this identical fix
      # was applied here too. The same derivation as Step 5 applies with
      # sum_term in place of quad_form and G in place of G_net. Skipped
      # entirely if use_horseshoe_layer = FALSE, same rationale as Step 6.
      # ---------------------------------------------------------------- #
      if (use_horseshoe_layer) {
        sum_term <- sum(beta_tilde[, d]^2 / lambda2[, d])
        tau_zeta <- rgamma(1, shape = gigg_a + G / 2, rate = gigg_b + sum_term / 2)
        zeta2[d] <- 1 / tau_zeta
      }

      # ---------------------------------------------------------------- #
      # Step 8: residual variance sigma2_d
      # ---------------------------------------------------------------- #
      if (use_group_layer) {
        mu_final <- alpha[grp_d, d] + rowSums(M_d * beta[grp_d, , d])
      } else {
        mu_final <- alpha[grp_d, d] + as.numeric(M_d %*% beta_tilde[, d])
      }
      ss_resid <- sum((y_d - mu_final)^2)
      sigma2[d] <- rinvgamma(1, shape = 0.01 + length(rows_d) / 2,
                              rate = 0.01 + ss_resid / 2)
    }

    if (verbose > 0 && t %% verbose == 0) {
      cat(sprintf("[gibbs_sampler] iteration %d / %d\n", t, n_iter))
    }

    if (t %in% keep_iters) {
      keep_pos <- keep_pos + 1L
      out_alpha[, , keep_pos] <- alpha
      out_beta_tilde[, , keep_pos] <- beta_tilde
      out_sigma2[, keep_pos] <- sigma2
      # Store NA (rather than the unused initialization value) for any
      # layer disabled via use_network_layer / use_horseshoe_layer, so a
      # saved ablation fit cannot be misread as having an active network
      # or horseshoe layer at some particular (meaningless) value -- the
      # config list above already records which layers were active, but
      # this makes the parameter arrays themselves self-consistent too.
      out_kappa2[, keep_pos] <- if (use_network_layer) kappa2 else NA
      out_lambda2[, , keep_pos] <- if (use_horseshoe_layer) lambda2 else NA
      out_zeta2[, keep_pos] <- if (use_horseshoe_layer) zeta2 else NA
    }
  }

  list(
    alpha = out_alpha,
    beta_tilde = out_beta_tilde,
    sigma2 = out_sigma2,
    kappa2 = out_kappa2,
    lambda2 = out_lambda2,
    zeta2 = out_zeta2,
    gene_names = gene_names,
    group_levels = levels(factor(group)),
    keep_iters = keep_iters,
    # Layer-usage metadata, so a saved fit is self-documenting about which
    # ablation configuration (if any) produced it -- relevant when
    # comparing multiple fits later via cross-validated predictive
    # likelihood (pathway.tex Phase 3).
    config = list(
      use_group_layer = use_group_layer,
      use_network_layer = use_network_layer,
      use_horseshoe_layer = use_horseshoe_layer,
      delta = delta,
      gigg_a = gigg_a,
      gigg_b = gigg_b,
      flat_prior_var = flat_prior_var
    )
  )
}
