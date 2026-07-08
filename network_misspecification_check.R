#!/usr/bin/env Rscript
#
# network_misspecification_check.R
#
# Network-misspecification robustness check, following Pham, Carvalho,
# Schaus & Kolaczyk (2018), Section 5.2: randomly reassign a fraction of
# genes' positions in the pathway-derived adjacency matrix W (simulating
# the pathway database having incorrectly assigned those genes' pathway
# membership), refit the model, and confirm predictive performance
# degrades gracefully as the corruption fraction increases, rather than
# collapsing catastrophically at some threshold. This is the Phase 2
# validation item specified in pathway.tex Section 3 ("Network-misspecification
# robustness check") and tracked as outstanding in the project progress log.
#
# Design notes:
#   - Corruption is implemented as a relabeling of W's row/column names
#     among a randomly chosen subset of genes (size = frac_corrupt * G),
#     NOT a structural edit to W's values. This means the network's
#     degree distribution and edge-weight structure are held fixed; only
#     the correspondence between gene identity and network position is
#     corrupted -- directly analogous to Pham et al.'s own procedure of
#     randomly reassigning genes to incorrect pathways while leaving the
#     pathway network's own structure untouched.
#   - Positive-definiteness of (D_W - delta*W) is verified at every
#     corruption level before fitting (relabeling does not change the
#     matrix's eigenvalues mathematically, but this is verified directly
#     rather than assumed, consistent with every other use of this check
#     elsewhere in the project).
#   - Performance is measured via held-out (out-of-sample) prediction:
#     a random 80% of cell lines are used to fit the model, the remaining
#     20% are predicted and scored via RMSE, at each corruption level.
#     This is a genuine generalization check (unlike the in-sample
#     comparison used as a quick sanity check for the ablation
#     comparators), since a network corrupted enough to actively mislead
#     the model should show worse OUT-OF-SAMPLE performance specifically.
#   - A modest drug subset and chain length are used by default (not the
#     full 295-drug, multi-thousand-iteration production settings), since
#     this check needs to be run at 5+ corruption levels x several
#     replicates each -- full production scale at every combination would
#     take many times longer than the single production run already
#     completed, for a check whose purpose is qualitative (does
#     performance degrade gracefully) rather than a result to report as
#     the paper's primary finding.

source("gibbs_sampler.R")

# --------------------------------------------------------------------- #
# Step 1: Load real data
# --------------------------------------------------------------------- #

mut <- read.csv("mutation_matrix.csv", row.names = 1, check.names = FALSE)
resp <- read.csv("response_matrix.csv")
W <- as.matrix(read.csv("gene_adjacency_W.csv", row.names = 1, check.names = FALSE))

model_ids <- colnames(mut)
M_full <- t(as.matrix(mut))
rownames(M_full) <- model_ids

cat(sprintf("Loaded: %d genes x %d models, network %d x %d\n",
            nrow(mut), ncol(mut), nrow(W), ncol(W)))

# --------------------------------------------------------------------- #
# Step 2: Corruption function
# --------------------------------------------------------------------- #

corrupt_network <- function(W, frac_corrupt, seed) {
  set.seed(seed)
  G <- nrow(W)
  n_corrupt <- round(frac_corrupt * G)
  if (n_corrupt < 2) return(W)  # need >=2 genes to swap labels meaningfully
  corrupt_idx <- sample(seq_len(G), n_corrupt)
  perm <- sample(corrupt_idx)
  new_names <- rownames(W)
  new_names[corrupt_idx] <- rownames(W)[perm]
  W_corrupt <- W
  rownames(W_corrupt) <- new_names
  colnames(W_corrupt) <- new_names
  W_corrupt
}

verify_pd <- function(W, delta = 0.8) {
  D_W <- diag(rowSums(W))
  eigvals <- eigen(D_W - delta * W, symmetric = TRUE, only.values = TRUE)$values
  min(eigvals) > 1e-8
}

# --------------------------------------------------------------------- #
# Step 3: Train/test split (held out at the CELL-LINE level, applied
# identically across all corruption levels and replicates so results are
# directly comparable)
# --------------------------------------------------------------------- #

set.seed(2024)
n_models <- length(model_ids)
test_frac <- 0.2
test_idx <- sample(seq_len(n_models), round(test_frac * n_models))
train_idx <- setdiff(seq_len(n_models), test_idx)
cat(sprintf("Train/test split: %d train, %d test cell lines\n",
            length(train_idx), length(test_idx)))

# --------------------------------------------------------------------- #
# Step 4: Fit-and-score function for one (corruption level, drug subset)
# combination
# --------------------------------------------------------------------- #

fit_and_score <- function(W_use, drug_ids, n_iter, burn_in, thin, seed) {
  resp_sub <- resp[resp$drug_id %in% drug_ids, ]
  y_full <- matrix(NA, nrow = n_models, ncol = length(drug_ids),
                    dimnames = list(model_ids, as.character(drug_ids)))
  for (i in seq_len(nrow(resp_sub))) {
    y_full[as.character(resp_sub$model_id[i]), as.character(resp_sub$drug_id[i])] <- resp_sub$y[i]
  }

  y_train <- y_full
  y_train[test_idx, ] <- NA  # mask test rows during fitting

  fit <- gibbs_sampler(
    y = y_train, M = M_full, group = rep(1, n_models), W = W_use,
    delta = 0.8, use_group_layer = FALSE,
    n_iter = n_iter, burn_in = burn_in, thin = thin, verbose = 0, seed = seed
  )

  beta_tilde_mean <- apply(fit$beta_tilde, c(1, 2), mean)
  alpha_mean <- apply(fit$alpha, c(1, 2), mean)

  sq_err <- c()
  for (d in seq_along(drug_ids)) {
    obs_test <- test_idx[!is.na(y_full[test_idx, d])]
    if (length(obs_test) == 0) next
    mu_pred <- alpha_mean[1, d] + as.numeric(M_full[obs_test, , drop = FALSE] %*% beta_tilde_mean[, d])
    sq_err <- c(sq_err, (y_full[obs_test, d] - mu_pred)^2)
  }

  # In addition to held-out RMSE (an AGGREGATE measure across all genes
  # and test cell lines), track TP53's own coefficient and kappa2's
  # posterior mean directly. This distinguishes two different things a
  # flat RMSE curve could mean: (a) the model is genuinely robust to
  # network corruption because the data term dominates aggregate
  # predictive fit (a real, reportable finding), vs (b) the check is
  # underpowered to detect anything because the network layer's
  # influence is too small to show up in EITHER condition at this drug
  # count / chain length. If TP53's coefficient or kappa2 visibly responds
  # to corruption even while RMSE stays flat, that supports (a); if
  # neither responds, that points to (b) and a larger drug count / longer
  # chain should be used before trusting the RMSE result.
  tp53_col <- which(colnames(M_full) == "TP53")
  tp53_beta <- if (length(tp53_col) == 1) beta_tilde_mean[tp53_col, 1] else NA
  kappa2_mean <- mean(fit$kappa2[1, ], na.rm = TRUE)

  list(rmse = sqrt(mean(sq_err)), tp53_beta = tp53_beta, kappa2_mean = kappa2_mean)
}

# --------------------------------------------------------------------- #
# Step 5: Run the check across corruption levels
# --------------------------------------------------------------------- #

corruption_levels <- c(0.00, 0.01, 0.05, 0.10, 0.20, 0.32)
n_replicates <- as.integer(Sys.getenv("MISSPEC_REPS", "3"))
test_drugs <- sort(unique(resp$drug_id))[1:10]  # modest drug subset; see header note
n_iter <- 800; burn_in <- 200; thin <- 2

cat(sprintf("\nRunning robustness check: %d corruption levels x %d replicates, %d drugs, n_iter=%d\n",
            length(corruption_levels), n_replicates, length(test_drugs), n_iter))

results <- data.frame(frac_corrupt = numeric(0), replicate = integer(0),
                       rmse = numeric(0), tp53_beta = numeric(0), kappa2_mean = numeric(0))

for (frac in corruption_levels) {
  for (rep_i in seq_len(n_replicates)) {
    W_use <- if (frac == 0) W else corrupt_network(W, frac, seed = 1000 + rep_i)
    pd_ok <- verify_pd(W_use)
    if (!pd_ok) {
      cat(sprintf("  WARNING: frac=%.2f rep=%d failed positive-definiteness check -- skipped\n", frac, rep_i))
      next
    }
    out <- fit_and_score(W_use, test_drugs, n_iter, burn_in, thin, seed = rep_i)
    cat(sprintf("  frac_corrupt=%.2f  rep=%d  held-out RMSE=%.4f  TP53 beta=%.4f  kappa2=%.5f\n",
                frac, rep_i, out$rmse, out$tp53_beta, out$kappa2_mean))
    results <- rbind(results, data.frame(frac_corrupt = frac, replicate = rep_i,
                                          rmse = out$rmse, tp53_beta = out$tp53_beta,
                                          kappa2_mean = out$kappa2_mean))
  }
}

# --------------------------------------------------------------------- #
# Step 6: Summarize and assess "graceful" vs "catastrophic" degradation
# --------------------------------------------------------------------- #

summary_df <- aggregate(rmse ~ frac_corrupt, data = results, FUN = function(x) c(mean = mean(x), sd = sd(x)))
summary_df <- do.call(data.frame, summary_df)
names(summary_df) <- c("frac_corrupt", "mean_rmse", "sd_rmse")
cat("\n=== Summary: mean held-out RMSE by corruption level ===\n")
print(summary_df)

cat("\n=== Gene-level sensitivity check (TP53 coefficient, kappa2) ===\n")
gene_sensitivity <- aggregate(cbind(tp53_beta, kappa2_mean) ~ frac_corrupt, data = results, FUN = mean)
print(gene_sensitivity)
cat("\nIf RMSE above is flat/non-monotonic across corruption levels, check whether\n")
cat("tp53_beta and kappa2_mean ALSO stay flat (suggesting the check is underpowered\n")
cat("to detect network corruption at this drug count/chain length) or whether they\n")
cat("visibly shift despite flat RMSE (suggesting the model is genuinely robust in\n")
cat("aggregate predictive fit despite real gene-level sensitivity to the network --\n")
cat("a meaningful, reportable distinction, not just noise).\n")

baseline_rmse <- summary_df$mean_rmse[summary_df$frac_corrupt == 0]
summary_df$pct_increase_vs_baseline <- 100 * (summary_df$mean_rmse - baseline_rmse) / baseline_rmse
cat("\n=== Percent RMSE increase relative to uncorrupted network ===\n")
print(summary_df[, c("frac_corrupt", "pct_increase_vs_baseline")])

# A simple, auditable "graceful degradation" check: RMSE should be
# monotonically non-decreasing (within noise) as corruption increases,
# and should NOT show a large jump (>50% relative increase) between any
# two adjacent corruption levels, which would indicate a discontinuous
# (catastrophic) failure rather than a smooth degradation.
cat("\n=== Graceful-degradation check ===\n")
max_jump <- max(diff(summary_df$mean_rmse) / summary_df$mean_rmse[-nrow(summary_df)]) * 100
cat(sprintf("Largest single-step relative RMSE increase between adjacent corruption levels: %.1f%%\n", max_jump))
if (max_jump > 50) {
  cat("FLAG: a jump exceeding 50% suggests a possible discontinuous (catastrophic) failure point -- inspect the corruption levels around this jump directly.\n")
} else {
  cat("No jump exceeds 50% -- consistent with graceful (smooth) degradation rather than catastrophic failure.\n")
}

write.csv(results, "network_misspecification_results.csv", row.names = FALSE)
write.csv(summary_df, "network_misspecification_summary.csv", row.names = FALSE)
cat("\nSaved network_misspecification_results.csv and network_misspecification_summary.csv\n")
