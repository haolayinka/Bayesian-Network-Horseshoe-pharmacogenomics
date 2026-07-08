#!/usr/bin/env Rscript

# =============================================================================
# Real Data Group-Layer Gibbs Sampler Runner with Targeted beta_group Storage
#
# Purpose:
#   Refit the revised group-layer model while storing posterior draws of
#   beta_group[k, g, d] for a targeted set of candidate genes. This is intended
#   for formal tissue/group-specific posterior summaries and credible intervals
#   in the manuscript group-layer effect section.
#
# Usage:
#   Rscript real_data_group_layer_vs_target_storage.R full_group_vs_target
#   Rscript real_data_group_layer_vs_target_storage.R nonetwork_group_vs_target
#   Rscript real_data_group_layer_vs_target_storage.R nohorseshoe_group_vs_target
#
# Recommended manuscript run:
#   Rscript real_data_group_layer_vs_target_storage.R full_group_vs_target
#
# Output files:
#   gibbs_fits_full_group_vs_target.rds
#   gibbs_fits_nonetwork_group_vs_target.rds
#   gibbs_fits_nohorseshoe_group_vs_target.rds
#
# Targeted storage added:
#   beta_group_target      [K x T x D x n_keep]
#   beta_group_target_mean [K x T x D]
#   target_gene_names      length T
#   target_gene_idx        indices in full gene panel
#
# Notes:
#   MIN_GROUP_SIZE is kept at 50 because the previous n_min=30 group-layer
#   results were less stable/less interpretable, while n_min=50 yielded K=6
#   larger groups and more stable group-level updates.
# =============================================================================

# -- 0. Configuration ----------------------------------------------------------

# File paths -- adjust as needed
MUTATION_FILE  <- "mutation_matrix.csv"
RESPONSE_FILE  <- "response_matrix.csv"
NETWORK_FILE   <- "gene_adjacency_W.csv"
GROUPS_FILE    <- "model_groups.csv"
OUTPUT_PREFIX  <- "gibbs_fits"

# Preprocessing
MIN_GROUP_SIZE <- 50        # Groups below this are merged into "Other"

# Model
DELTA          <- 0.8       # Network smoothing, fixed
RIDGE_EPS      <- 1e-4      # Ridge added to A_k

# GIGG hyperparameters for kappa2 and zeta2 precisions
GIGG_A         <- 1.0
GIGG_B         <- 1.0

# eta2 hyperparameters: proper prior replacing IG(0.01, 0.01)
ETA2_A         <- 2.0       # IG shape: mode = b/(a+1) = 0.5/3 = 0.167
ETA2_B         <- 0.5       # IG rate:  mean = b/(a-1) = 0.500

# sigma2 hyperparameters
SIGMA2_A       <- 0.01
SIGMA2_B       <- 0.01

# alpha prior variance
ALPHA_VAR      <- 100

# MCMC
N_ITER         <- 5000
BURN_IN        <- 1000
THIN           <- 2
N_CHAINS       <- 2

# Target genes for full group-specific posterior storage.
# These include the main no-group full-model findings and biologically coherent
# group-layer signals from the exploratory VS scan.
TARGET_GENES <- unique(c(
  "BRAF", "EGFR", "TP53", "ASXL1", "STAG2", "PTEN", "KMT2D", "EZH2",
  "PBRM1", "RB1", "KRAS", "NF1", "TET2", "NRAS",
  "APC", "PCDH17", "PRKD2", "ZNF626", "EBF1", "CACNA1D",
  "STAT3", "ATM", "STK11", "ALK", "BAP1"
))

cat("Configuration:\n")
cat(sprintf("  MIN_GROUP_SIZE: %d\n", MIN_GROUP_SIZE))
cat(sprintf("  DELTA:          %.2f\n", DELTA))
cat(sprintf("  RIDGE_EPS:      %g\n", RIDGE_EPS))
cat(sprintf("  ETA2 prior:     IG(%.1f, %.1f) [mode=%.3f, mean=%.3f]\n",
            ETA2_A, ETA2_B,
            ETA2_B / (ETA2_A + 1),
            ETA2_B / (ETA2_A - 1)))
cat(sprintf("  GIGG prior:     Gamma(%.1f, %.1f) on precisions\n", GIGG_A, GIGG_B))
cat(sprintf("  MCMC:           %d iter / %d burn-in / thin=%d / %d chains\n",
            N_ITER, BURN_IN, THIN, N_CHAINS))
cat(sprintf("  Target genes requested: %d\n\n", length(TARGET_GENES)))

# -- 1. Load data --------------------------------------------------------------

cat("Loading data...\n")
mut       <- read.csv(MUTATION_FILE, row.names = 1, check.names = FALSE)
resp      <- read.csv(RESPONSE_FILE)
W_raw     <- read.csv(NETWORK_FILE, row.names = 1, check.names = FALSE)
groups_df <- read.csv(GROUPS_FILE)

gene_names     <- rownames(mut)
model_ids_all  <- colnames(mut)
G              <- length(gene_names)
all_drugs      <- sort(unique(resp$drug_id))
D              <- length(all_drugs)

cat(sprintf("  Genes: %d | Drugs: %d | Cell lines (raw): %d\n",
            G, D, length(model_ids_all)))

# Target gene indices in full panel
target_idx <- match(TARGET_GENES, gene_names)
missing_targets <- TARGET_GENES[is.na(target_idx)]
target_idx <- target_idx[!is.na(target_idx)]
target_gene_names <- gene_names[target_idx]
T_TARGET <- length(target_idx)

cat(sprintf("  Target genes matched in panel: %d / %d\n", T_TARGET, length(TARGET_GENES)))
if (length(missing_targets) > 0) {
  cat("  Target genes not found and skipped:\n")
  print(missing_targets)
}
if (T_TARGET == 0) stop("No target genes were found in gene_names.")
cat("  Target genes stored:\n")
print(target_gene_names)

# -- 2. Tissue group assignment and merging -----------------------------------

cat(sprintf("\nMerging tissue groups with n < %d into 'Other'...\n", MIN_GROUP_SIZE))

group_raw <- groups_df$tissue[match(model_ids_all, groups_df$model_id)]

tissue_counts <- sort(table(group_raw[!is.na(group_raw)]))
cat("Group sizes before merging:\n")
print(tissue_counts)

group_merged <- ifelse(
  is.na(group_raw),
  "Other",
  ifelse(tissue_counts[group_raw] < MIN_GROUP_SIZE, "Other", group_raw)
)

valid_cells <- !is.na(group_raw)
cat(sprintf("\nCell lines with tissue labels: %d / %d\n",
            sum(valid_cells), length(valid_cells)))

model_ids <- model_ids_all[valid_cells]
group_vec <- group_merged[valid_cells]
mut_use   <- mut[, valid_cells, drop = FALSE]

tissue_levels <- sort(unique(group_vec))
K             <- length(tissue_levels)
group_idx     <- match(group_vec, tissue_levels)

cat(sprintf("\nAfter merging: K = %d tissue groups\n", K))
cat("Group sizes after merging:\n")
print(sort(table(group_vec)))

N <- length(model_ids)

# -- 3. Network matrix ---------------------------------------------------------

cat("\nPreparing network matrix...\n")
net_genes <- intersect(gene_names, rownames(W_raw))
W_net     <- as.matrix(W_raw[net_genes, net_genes])
net_idx   <- match(net_genes, gene_names)
G_net     <- length(net_genes)

D_W <- diag(rowSums(W_net))
Q   <- D_W - DELTA * W_net

eigs <- eigen(Q, only.values = TRUE)$values
if (any(eigs <= 0)) {
  stop(sprintf("Q is not positive-definite at delta=%.1f. Minimum eigenvalue: %.6f",
               DELTA, min(eigs)))
}
cat(sprintf("  Network: %d genes | delta=%.1f | min eigenvalue=%.4f (PD OK)\n",
            G_net, DELTA, min(eigs)))
Prec_net_raw <- Q

# -- 4. Gibbs sampler function -------------------------------------------------

run_gibbs <- function(chain_id,
                      seed_id = chain_id,
                      use_network = TRUE,
                      use_horseshoe = TRUE) {

  cat(sprintf("\n=== Chain %d | seed=%d | network=%s | horseshoe=%s ===\n",
              chain_id, seed_id, use_network, use_horseshoe))
  set.seed(seed_id)

  keep_every <- THIN
  n_keep     <- floor((N_ITER - BURN_IN) / keep_every)
  keep_iters <- seq(BURN_IN + keep_every, N_ITER, by = keep_every)

  # -- Storage ----------------------------------------------------------------
  beta_tilde_store <- array(0, dim = c(G, D, n_keep))

  # Full beta_group posterior mean only, as before.
  beta_group_sum <- array(0, dim = c(K, G, D))

  # NEW: targeted full posterior draws of beta_group for candidate genes.
  beta_group_target_store <- array(0, dim = c(K, T_TARGET, D, n_keep))
  beta_group_target_sum   <- array(0, dim = c(K, T_TARGET, D))

  keep_pos <- 0L

  alpha_store  <- array(0, dim = c(K, D, n_keep))
  sigma2_store <- matrix(0, D, n_keep)
  kappa2_store <- matrix(0, D, n_keep)
  zeta2_store  <- matrix(0, D, n_keep)

  # -- Initialize parameters --------------------------------------------------
  beta_tilde <- matrix(0, G, D)
  beta_group <- array(0, dim = c(K, G, D))
  alpha      <- matrix(0, K, D)
  sigma2     <- rep(1, D)
  kappa2     <- rep(0.5, D)
  zeta2      <- rep(0.5, D)
  eta2       <- matrix(ETA2_B / (ETA2_A + 1), G, D)
  lambda2    <- matrix(1, G, D)
  nu         <- matrix(1, G, D)

  # -- Precompute per-drug data structures -----------------------------------
  cat("  Precomputing per-drug structures...\n")
  drug_data <- vector("list", D)
  for (d in seq_len(D)) {
    drug_id <- all_drugs[d]
    resp_d  <- resp[resp$drug_id == drug_id, ]
    obs_idx <- match(resp_d$model_id, model_ids)
    valid   <- !is.na(obs_idx)
    obs_idx <- obs_idx[valid]
    y_d     <- resp_d$y[valid]
    if (length(obs_idx) == 0) next

    M_d   <- t(as.matrix(mut_use[, obs_idx, drop = FALSE]))
    grp_d <- group_idx[obs_idx]

    grp_list <- lapply(seq_len(K), function(k) which(grp_d == k))

    drug_data[[d]] <- list(
      obs_idx  = obs_idx,
      y        = y_d,
      M        = M_d,
      grp      = grp_d,
      grp_list = grp_list,
      n        = length(y_d)
    )
  }
  cat("  Starting MCMC...\n")

  pb_step <- max(1, floor(N_ITER / 20))
  for (iter in seq_len(N_ITER)) {

    if (iter %% pb_step == 0) {
      cat(sprintf("    iter %4d / %d\n", iter, N_ITER))
      flush.console()
    }

    for (d in seq_len(D)) {
      dd <- drug_data[[d]]
      if (is.null(dd)) next

      y_d      <- dd$y
      M_d      <- dd$M
      grp_d    <- dd$grp
      grp_list <- dd$grp_list
      n_d      <- dd$n

      # Step 1: alpha[k, d]
      for (k in seq_len(K)) {
        rows_k <- grp_list[[k]]
        if (length(rows_k) == 0) {
          alpha[k, d] <- rnorm(1, 0, sqrt(ALPHA_VAR))
          next
        }
        mu_mut_k <- rowSums(
          M_d[rows_k, , drop = FALSE] *
            matrix(beta_group[k, , d], nrow = length(rows_k), ncol = G, byrow = TRUE)
        )
        resid_k <- y_d[rows_k] - mu_mut_k
        prec_k  <- length(rows_k) / sigma2[d] + 1 / ALPHA_VAR
        mean_k  <- sum(resid_k) / sigma2[d] / prec_k
        alpha[k, d] <- rnorm(1, mean_k, 1 / sqrt(prec_k))
      }

      # Step 2: beta_group[k, , d]
      for (k in seq_len(K)) {
        rows_k <- grp_list[[k]]
        if (length(rows_k) == 0) {
          beta_group[k, , d] <- rnorm(G, beta_tilde[, d], sqrt(eta2[, d]))
          next
        }
        y_k <- y_d[rows_k] - alpha[k, d]
        M_k <- M_d[rows_k, , drop = FALSE]

        MtM_k <- crossprod(M_k) / sigma2[d]
        A_k   <- MtM_k + diag(1 / eta2[, d]) + RIDGE_EPS * diag(G)
        b_k   <- as.numeric(crossprod(M_k, y_k)) / sigma2[d] +
                 beta_tilde[, d] / eta2[, d]

        R_k <- tryCatch(chol(A_k), error = function(e) {
          chol(A_k + 0.1 * diag(G))
        })
        mean_k <- backsolve(R_k, forwardsolve(t(R_k), b_k))
        draw_k <- backsolve(R_k, rnorm(G))
        beta_group[k, , d] <- mean_k + draw_k
      }

      # Step 3: eta2[g, d]
      for (g in seq_len(G)) {
        ss_g    <- sum((beta_group[, g, d] - beta_tilde[g, d])^2)
        shape_g <- ETA2_A + K / 2
        rate_g  <- ETA2_B + ss_g / 2
        eta2[g, d] <- 1 / rgamma(1, shape = shape_g, rate = rate_g)
      }

      # Step 4: beta_tilde[, d]
      data_prec_diag <- rep(0, G)
      data_b         <- rep(0, G)
      for (k in seq_len(K)) {
        data_prec_diag <- data_prec_diag + 1 / eta2[, d]
        data_b         <- data_b + beta_group[k, , d] / eta2[, d]
      }

      net_prec_full <- matrix(0, G, G)
      if (use_network) {
        net_prec_full[net_idx, net_idx] <- Prec_net_raw / kappa2[d]
      }

      if (use_horseshoe) {
        hs_prec_diag <- 1 / (lambda2[, d] * zeta2[d])
      } else {
        hs_prec_diag <- rep(1 / 100, G)
      }

      B_d <- net_prec_full + diag(data_prec_diag + hs_prec_diag)
      b_d <- data_b

      R_d     <- chol(B_d)
      mean_bt <- backsolve(R_d, forwardsolve(t(R_d), b_d))
      draw_bt <- backsolve(R_d, rnorm(G))
      beta_tilde[, d] <- mean_bt + draw_bt

      # Step 5: kappa2[d]
      if (use_network) {
        bt_net <- beta_tilde[net_idx, d]
        Q_d    <- as.numeric(t(bt_net) %*% Prec_net_raw %*% bt_net)
        tau_k  <- rgamma(1, shape = GIGG_A + G_net / 2,
                         rate = GIGG_B + Q_d / 2)
        kappa2[d] <- 1 / tau_k
      }

      # Step 6: lambda2[g, d] and nu[g, d]
      if (use_horseshoe) {
        for (g in seq_len(G)) {
          nu[g, d] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda2[g, d])
          rate_lam <- 1 / nu[g, d] + beta_tilde[g, d]^2 / (2 * zeta2[d])
          lambda2[g, d] <- 1 / rgamma(1, shape = 1, rate = rate_lam)
        }
      }

      # Step 7: zeta2[d]
      if (use_horseshoe) {
        S_d   <- sum(beta_tilde[, d]^2 / lambda2[, d])
        tau_z <- rgamma(1, shape = GIGG_A + G / 2,
                        rate = GIGG_B + S_d / 2)
        zeta2[d] <- 1 / tau_z
      }

      # Step 8: sigma2[d]
      beta_group_rows <- beta_group[grp_d, , d]
      if (is.null(dim(beta_group_rows))) {
        beta_group_rows <- matrix(beta_group_rows, nrow = nrow(M_d), ncol = G, byrow = FALSE)
      }

      mu_full <- alpha[grp_d, d] + rowSums(M_d * beta_group_rows)

      ss_resid  <- sum((y_d - mu_full)^2)
      sigma2[d] <- 1 / rgamma(1,
                              shape = SIGMA2_A + n_d / 2,
                              rate  = SIGMA2_B + ss_resid / 2)
    }

    # Store if past burn-in
    if (iter > BURN_IN && (iter - BURN_IN) %% keep_every == 0) {
      keep_pos <- keep_pos + 1L
      s        <- keep_pos

      beta_tilde_store[, , s] <- beta_tilde
      alpha_store[, , s]      <- alpha
      sigma2_store[, s]       <- sigma2
      kappa2_store[, s]       <- kappa2
      zeta2_store[, s]        <- zeta2

      beta_group_sum <- beta_group_sum + beta_group

      # NEW targeted group-specific posterior draws and running mean.
      beta_group_target_store[, , , s] <- beta_group[, target_idx, , drop = FALSE]
      beta_group_target_sum <- beta_group_target_sum + beta_group[, target_idx, , drop = FALSE]
    }
  }

  list(
    beta_tilde             = beta_tilde_store,
    beta_group_mean        = beta_group_sum / keep_pos,
    beta_group_target      = beta_group_target_store,
    beta_group_target_mean = beta_group_target_sum / keep_pos,
    target_gene_names      = target_gene_names,
    target_gene_idx        = target_idx,
    alpha                  = alpha_store,
    sigma2                 = sigma2_store,
    kappa2                 = kappa2_store,
    zeta2                  = zeta2_store,
    gene_names             = gene_names,
    group_levels           = tissue_levels,
    keep_iters             = keep_iters,
    config = list(
      use_network        = use_network,
      use_horseshoe      = use_horseshoe,
      delta              = DELTA,
      ridge_eps          = RIDGE_EPS,
      eta2_a             = ETA2_A,
      eta2_b             = ETA2_B,
      gigg_a             = GIGG_A,
      gigg_b             = GIGG_B,
      min_group_size     = MIN_GROUP_SIZE,
      n_iter             = N_ITER,
      burn_in            = BURN_IN,
      thin               = THIN,
      K                  = K,
      tissue_levels      = tissue_levels,
      target_gene_names  = target_gene_names,
      target_gene_idx    = target_idx,
      target_storage_dim = dim(beta_group_target_store)
    )
  )
}

# -- 5. Fit one model configuration from command-line argument -----------------

args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "help")

cat("\nRequested mode:", mode, "\n")

model_configs <- list(
  full_group_vs_target = list(
    name = "full_group_vs_target",
    network = TRUE,
    horseshoe = TRUE,
    seed_offset = 10000
  ),
  nonetwork_group_vs_target = list(
    name = "nonetwork_group_vs_target",
    network = FALSE,
    horseshoe = TRUE,
    seed_offset = 20000
  ),
  nohorseshoe_group_vs_target = list(
    name = "nohorseshoe_group_vs_target",
    network = TRUE,
    horseshoe = FALSE,
    seed_offset = 30000
  )
)

if (mode == "help" || !mode %in% names(model_configs)) {
  cat("\nUsage:\n")
  cat("  Rscript real_data_group_layer_vs_target_storage.R full_group_vs_target\n")
  cat("  Rscript real_data_group_layer_vs_target_storage.R nonetwork_group_vs_target\n")
  cat("  Rscript real_data_group_layer_vs_target_storage.R nohorseshoe_group_vs_target\n\n")
  cat("Recommended manuscript run:\n")
  cat("  Rscript real_data_group_layer_vs_target_storage.R full_group_vs_target\n\n")
  quit(save = "no", status = 0)
}

cfg <- model_configs[[mode]]

cat("\nRunning model:\n")
print(cfg)

fits <- vector("list", N_CHAINS)

for (ch in seq_len(N_CHAINS)) {
  fits[[ch]] <- run_gibbs(
    chain_id = ch,
    seed_id = cfg$seed_offset + ch,
    use_network = cfg$network,
    use_horseshoe = cfg$horseshoe
  )

  out_file_partial <- sprintf("%s_%s_partial.rds", OUTPUT_PREFIX, cfg$name)
  saveRDS(fits, out_file_partial)
  cat(sprintf("\nSaved partial file after chain %d: %s\n", ch, out_file_partial))
}

out_file <- sprintf("%s_%s.rds", OUTPUT_PREFIX, cfg$name)
saveRDS(fits, out_file)

cat(sprintf("\nSaved final file: %s\n", out_file))

cat("\nFit summary:\n")
cat("Number of chains:", length(fits), "\n")
cat("NULL chains:\n")
print(sapply(fits, is.null))

cat("beta_tilde dimensions:\n")
print(dim(fits[[1]]$beta_tilde))
print(dim(fits[[2]]$beta_tilde))

cat("beta_group_mean dimensions:\n")
print(dim(fits[[1]]$beta_group_mean))
print(dim(fits[[2]]$beta_group_mean))

cat("beta_group_target dimensions:\n")
print(dim(fits[[1]]$beta_group_target))
print(dim(fits[[2]]$beta_group_target))

cat("Target genes stored:\n")
print(fits[[1]]$target_gene_names)

cat("\n=== Fit complete ===\n")
