#!/usr/bin/env Rscript
# ============================================================
# scan_tissue_effects_vs_target_posterior.R
#
# Uses the NEW targeted-storage group-layer run:
#   gibbs_fits_full_group_vs_target.rds
#
# This script computes formal posterior summaries and credible intervals
# from stored posterior draws:
#   beta_group_target[k, target_gene, d, sample]
#
# It also computes tissue/group deviations from the shared effect:
#   beta_group_target[k, target_gene, d, sample]
#     - beta_tilde[target_gene_index, d, sample]
#
# This replaces the older exploratory scan that used beta_group_mean and
# between-chain SD. The intervals here are true posterior intervals based
# on retained MCMC draws for the target genes.
#
# Output CSV files:
#   tissue_target_ci_scan_full_vs.csv
#   tissue_target_beta_80_included_vs.csv
#   tissue_target_beta_90_included_vs.csv
#   tissue_target_beta_95_included_vs.csv
#   tissue_target_delta_80_included_vs.csv
#   tissue_target_delta_90_included_vs.csv
#   tissue_target_delta_95_included_vs.csv
#   tissue_target_gene_summary_vs.csv
#   tissue_target_ezh2_kmt2d_vs.csv
#   tissue_target_braf_egfr_vs.csv
#   tissue_target_top_beta_effects_vs.csv
#   tissue_target_top_deviations_vs.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

cat("=== Formal Tissue/Group-Specific Posterior Scan: Target Genes ===\n\n")

# -- 0. File paths -------------------------------------------------------------
FIT_FILE      <- "gibbs_fits_full_group_vs_target.rds"
RESPONSE_FILE <- "response_matrix.csv"
GROUPS_FILE   <- "model_groups.csv"
MUTATION_FILE <- "mutation_matrix.csv"
OUT_PREFIX    <- "tissue_target"

# -- 1. Load fit ---------------------------------------------------------------
cat(sprintf("Loading %s...\n", FIT_FILE))
if (!file.exists(FIT_FILE)) stop("Cannot find fit file: ", FIT_FILE)

fits <- readRDS(FIT_FILE)
n_chains <- length(fits)
if (n_chains < 1) stop("No chains found in fit object.")
if (any(sapply(fits, is.null))) stop("At least one chain is NULL.")

# Required objects
required_names <- c("beta_group_target", "beta_tilde", "target_gene_names",
                    "target_gene_idx", "group_levels", "gene_names", "config")
missing_names <- setdiff(required_names, names(fits[[1]]))
if (length(missing_names) > 0) {
  stop("Fit object is missing required entries: ", paste(missing_names, collapse = ", "))
}

TISSUE_LEVELS <- fits[[1]]$group_levels
K <- length(TISSUE_LEVELS)

target_gene_names <- fits[[1]]$target_gene_names
target_gene_idx <- fits[[1]]$target_gene_idx
T_TARGET <- length(target_gene_names)

gene_names <- fits[[1]]$gene_names
G <- length(gene_names)
D <- dim(fits[[1]]$beta_tilde)[2]
n_keep <- dim(fits[[1]]$beta_tilde)[3]

cat(sprintf("  Chains: %d | K=%d | target genes=%d | G=%d | D=%d | samples/chain=%d\n",
            n_chains, K, T_TARGET, G, D, n_keep))
cat("  Group levels:\n")
print(TISSUE_LEVELS)
cat("  Target genes:\n")
print(target_gene_names)

# Dimension checks for every chain
for (ch in seq_len(n_chains)) {
  d_bg <- dim(fits[[ch]]$beta_group_target)
  d_bt <- dim(fits[[ch]]$beta_tilde)
  cat(sprintf("  Chain %d beta_group_target dim: [%s]\n", ch, paste(d_bg, collapse = " x ")))
  cat(sprintf("  Chain %d beta_tilde dim:       [%s]\n", ch, paste(d_bt, collapse = " x ")))
  stopifnot(d_bg[1] == K, d_bg[2] == T_TARGET, d_bg[3] == D)
  stopifnot(d_bt[1] == G, d_bt[2] == D)
}

# -- 2. Load drug names --------------------------------------------------------
cat("\nLoading drug names...\n")
if (!file.exists(RESPONSE_FILE)) stop("Cannot find response file: ", RESPONSE_FILE)
resp <- read.csv(RESPONSE_FILE)
if (!"drug_id" %in% names(resp)) stop("response_matrix.csv must contain drug_id.")
if (!"drug_name" %in% names(resp)) {
  warning("response_matrix.csv does not contain drug_name; using drug_id as drug_name.")
  resp$drug_name <- as.character(resp$drug_id)
}

drug_lookup <- unique(resp[, c("drug_id", "drug_name")])
all_drugs <- sort(unique(resp$drug_id))
stopifnot(length(all_drugs) == D)

get_drug_name <- function(id) {
  out <- drug_lookup$drug_name[match(id, drug_lookup$drug_id)]
  ifelse(is.na(out), as.character(id), as.character(out))
}

# -- 3. Reconstruct group sizes for annotation --------------------------------
cat("\nReconstructing group sizes from model_groups.csv...\n")
n_k <- rep(NA_integer_, K)
names(n_k) <- TISSUE_LEVELS

if (file.exists(GROUPS_FILE) && file.exists(MUTATION_FILE)) {
  groups_df <- read.csv(GROUPS_FILE)
  mut <- read.csv(MUTATION_FILE, row.names = 1, check.names = FALSE)
  model_ids_all <- colnames(mut)

  MIN_GROUP_SIZE <- fits[[1]]$config$min_group_size
  if (is.null(MIN_GROUP_SIZE)) MIN_GROUP_SIZE <- 50L
  MIN_GROUP_SIZE <- as.integer(MIN_GROUP_SIZE)

  group_raw <- groups_df$tissue[match(model_ids_all, groups_df$model_id)]
  tc <- sort(table(group_raw[!is.na(group_raw)]))
  group_merged <- ifelse(
    is.na(group_raw),
    "Other",
    ifelse(tc[group_raw] < MIN_GROUP_SIZE, "Other", group_raw)
  )
  valid_cells <- !is.na(group_raw)
  group_vec <- group_merged[valid_cells]
  n_k <- as.integer(table(factor(group_vec, levels = TISSUE_LEVELS)))
  names(n_k) <- TISSUE_LEVELS
  cat(sprintf("  MIN_GROUP_SIZE: %d\n", MIN_GROUP_SIZE))
  cat("  Group sizes aligned to fitted group_levels:\n")
  print(n_k)
} else {
  warning("model_groups.csv or mutation_matrix.csv not found; n_k will be NA.")
}

# -- 4. Helper functions -------------------------------------------------------
qfun <- function(x, probs) as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))

summarise_draws <- function(x) {
  qs <- qfun(x, c(0.025, 0.05, 0.10, 0.90, 0.95, 0.975))
  p_pos <- mean(x > 0, na.rm = TRUE)
  p_neg <- mean(x < 0, na.rm = TRUE)
  c(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = qs[1],
    q05 = qs[2],
    q10 = qs[3],
    q90 = qs[4],
    q95 = qs[5],
    q975 = qs[6],
    p_pos = p_pos,
    p_neg = p_neg,
    sign_prob = max(p_pos, p_neg)
  )
}

# -- 5. Scan K x target genes x D ---------------------------------------------
cat("\n=== Scanning K x target genes x D combinations using posterior draws ===\n")
cat(sprintf("Total combinations: %d x %d x %d = %d\n\n",
            K, T_TARGET, D, K * T_TARGET * D))

rows <- vector("list", K * T_TARGET * D)
idx <- 1L

for (k in seq_len(K)) {
  cat(sprintf("  Group %d/%d: %s\n", k, K, TISSUE_LEVELS[k]))
  flush.console()

  for (tg in seq_len(T_TARGET)) {
    g_full <- target_gene_idx[tg]
    g_name <- target_gene_names[tg]

    for (d in seq_len(D)) {
      # Combine posterior draws across chains for group-specific effect
      beta_draws <- unlist(lapply(seq_len(n_chains), function(ch) {
        as.numeric(fits[[ch]]$beta_group_target[k, tg, d, ])
      }), use.names = FALSE)

      # Combine posterior draws across chains for shared effect, same gene/drug
      shared_draws <- unlist(lapply(seq_len(n_chains), function(ch) {
        as.numeric(fits[[ch]]$beta_tilde[g_full, d, ])
      }), use.names = FALSE)

      # Tissue/group-specific deviation from shared effect
      delta_draws <- beta_draws - shared_draws

      beta_sum <- summarise_draws(beta_draws)
      delta_sum <- summarise_draws(delta_draws)

      beta_sig80 <- (beta_sum["q10"] > 0) | (beta_sum["q90"] < 0)
      beta_sig90 <- (beta_sum["q05"] > 0) | (beta_sum["q95"] < 0)
      beta_sig95 <- (beta_sum["q025"] > 0) | (beta_sum["q975"] < 0)

      delta_sig80 <- (delta_sum["q10"] > 0) | (delta_sum["q90"] < 0)
      delta_sig90 <- (delta_sum["q05"] > 0) | (delta_sum["q95"] < 0)
      delta_sig95 <- (delta_sum["q025"] > 0) | (delta_sum["q975"] < 0)

      rows[[idx]] <- data.frame(
        tissue = TISSUE_LEVELS[k],
        n_k = n_k[k],
        gene = g_name,
        target_gene_index = tg,
        full_gene_index = g_full,
        drug_id = all_drugs[d],
        drug_name = get_drug_name(all_drugs[d]),

        beta_mean = beta_sum["mean"],
        beta_sd = beta_sum["sd"],
        beta_ci80_lo = beta_sum["q10"],
        beta_ci80_hi = beta_sum["q90"],
        beta_ci90_lo = beta_sum["q05"],
        beta_ci90_hi = beta_sum["q95"],
        beta_ci95_lo = beta_sum["q025"],
        beta_ci95_hi = beta_sum["q975"],
        beta_p_pos = beta_sum["p_pos"],
        beta_p_neg = beta_sum["p_neg"],
        beta_sign_prob = beta_sum["sign_prob"],
        beta_sig80 = beta_sig80,
        beta_sig90 = beta_sig90,
        beta_sig95 = beta_sig95,
        beta_direction = ifelse(beta_sum["mean"] > 0, "resistance", "sensitivity"),
        beta_abs_mean = abs(beta_sum["mean"]),

        delta_mean = delta_sum["mean"],
        delta_sd = delta_sum["sd"],
        delta_ci80_lo = delta_sum["q10"],
        delta_ci80_hi = delta_sum["q90"],
        delta_ci90_lo = delta_sum["q05"],
        delta_ci90_hi = delta_sum["q95"],
        delta_ci95_lo = delta_sum["q025"],
        delta_ci95_hi = delta_sum["q975"],
        delta_p_pos = delta_sum["p_pos"],
        delta_p_neg = delta_sum["p_neg"],
        delta_sign_prob = delta_sum["sign_prob"],
        delta_sig80 = delta_sig80,
        delta_sig90 = delta_sig90,
        delta_sig95 = delta_sig95,
        delta_direction = ifelse(delta_sum["mean"] > 0, "above_shared", "below_shared"),
        delta_abs_mean = abs(delta_sum["mean"]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
}

scan_df <- bind_rows(rows)
scan_df <- scan_df %>% arrange(desc(beta_abs_mean))

# -- 6. Flagged sets -----------------------------------------------------------
beta80  <- scan_df %>% filter(beta_sig80)
beta90  <- scan_df %>% filter(beta_sig90)
beta95  <- scan_df %>% filter(beta_sig95)
delta80 <- scan_df %>% filter(delta_sig80)
delta90 <- scan_df %>% filter(delta_sig90)
delta95 <- scan_df %>% filter(delta_sig95)

cat("\n=== Group-specific beta effect CI flags ===\n")
cat(sprintf("Beta flagged at 80%% CI: %d / %d (%.3f%%)\n",
            nrow(beta80), nrow(scan_df), 100 * nrow(beta80) / nrow(scan_df)))
cat(sprintf("Beta flagged at 90%% CI: %d / %d (%.3f%%)\n",
            nrow(beta90), nrow(scan_df), 100 * nrow(beta90) / nrow(scan_df)))
cat(sprintf("Beta flagged at 95%% CI: %d / %d (%.3f%%)\n",
            nrow(beta95), nrow(scan_df), 100 * nrow(beta95) / nrow(scan_df)))

cat("\n=== Tissue deviation from shared effect CI flags ===\n")
cat(sprintf("Delta flagged at 80%% CI: %d / %d (%.3f%%)\n",
            nrow(delta80), nrow(scan_df), 100 * nrow(delta80) / nrow(scan_df)))
cat(sprintf("Delta flagged at 90%% CI: %d / %d (%.3f%%)\n",
            nrow(delta90), nrow(scan_df), 100 * nrow(delta90) / nrow(scan_df)))
cat(sprintf("Delta flagged at 95%% CI: %d / %d (%.3f%%)\n",
            nrow(delta95), nrow(scan_df), 100 * nrow(delta95) / nrow(scan_df)))

# -- 7. Summaries --------------------------------------------------------------
cat("\n=== Top target genes by number of 95% beta-CI flagged tissue-drug pairs ===\n")
if (nrow(beta95) > 0) print(head(sort(table(beta95$gene), decreasing = TRUE), 25)) else cat("No beta 95% CI flags.\n")

cat("\n=== Top target genes by number of 95% delta-CI flagged tissue-drug deviations ===\n")
if (nrow(delta95) > 0) print(head(sort(table(delta95$gene), decreasing = TRUE), 25)) else cat("No delta 95% CI flags.\n")

cat("\n=== Top tissues by 95% beta-CI flagged count ===\n")
if (nrow(beta95) > 0) print(sort(table(beta95$tissue), decreasing = TRUE)) else cat("No beta 95% CI flags.\n")

cat("\n=== Top tissues by 95% delta-CI flagged count ===\n")
if (nrow(delta95) > 0) print(sort(table(delta95$tissue), decreasing = TRUE)) else cat("No delta 95% CI flags.\n")

gene_summary <- scan_df %>%
  group_by(gene) %>%
  summarise(
    n_beta80 = sum(beta_sig80),
    n_beta90 = sum(beta_sig90),
    n_beta95 = sum(beta_sig95),
    n_delta80 = sum(delta_sig80),
    n_delta90 = sum(delta_sig90),
    n_delta95 = sum(delta_sig95),
    mean_abs_beta = round(mean(beta_abs_mean), 4),
    max_abs_beta = round(max(beta_abs_mean), 4),
    mean_abs_delta = round(mean(delta_abs_mean), 4),
    max_abs_delta = round(max(delta_abs_mean), 4),
    top_beta_tissue = tissue[which.max(beta_abs_mean)],
    top_beta_drug = drug_name[which.max(beta_abs_mean)],
    top_beta_mean = beta_mean[which.max(beta_abs_mean)],
    top_delta_tissue = tissue[which.max(delta_abs_mean)],
    top_delta_drug = drug_name[which.max(delta_abs_mean)],
    top_delta_mean = delta_mean[which.max(delta_abs_mean)],
    .groups = "drop"
  ) %>%
  arrange(desc(n_beta95), desc(max_abs_beta))

cat("\n=== Gene summary ===\n")
print(gene_summary, n = Inf)

# Focused outputs
focus_genes_1 <- c("EZH2", "KMT2D")
focus_genes_2 <- c("BRAF", "EGFR", "TP53", "ASXL1", "NRAS", "STAG2")

ezh2_kmt2d <- scan_df %>%
  filter(gene %in% focus_genes_1) %>%
  arrange(gene, tissue, desc(beta_abs_mean))

braf_egfr <- scan_df %>%
  filter(gene %in% focus_genes_2) %>%
  arrange(gene, tissue, desc(beta_abs_mean))

top_beta_effects <- scan_df %>%
  arrange(desc(beta_abs_mean)) %>%
  select(tissue, n_k, gene, drug_id, drug_name,
         beta_mean, beta_sd, beta_ci80_lo, beta_ci80_hi,
         beta_ci95_lo, beta_ci95_hi, beta_sign_prob,
         beta_sig80, beta_sig90, beta_sig95, beta_direction,
         delta_mean, delta_ci95_lo, delta_ci95_hi, delta_sign_prob,
         delta_sig95) %>%
  head(50)

top_deviations <- scan_df %>%
  arrange(desc(delta_abs_mean)) %>%
  select(tissue, n_k, gene, drug_id, drug_name,
         delta_mean, delta_sd, delta_ci80_lo, delta_ci80_hi,
         delta_ci95_lo, delta_ci95_hi, delta_sign_prob,
         delta_sig80, delta_sig90, delta_sig95, delta_direction,
         beta_mean, beta_ci95_lo, beta_ci95_hi, beta_sign_prob,
         beta_sig95) %>%
  head(50)

cat("\n=== Top 20 group-specific beta effects by |posterior mean| ===\n")
print(head(top_beta_effects, 20), row.names = FALSE)

cat("\n=== Top 20 tissue deviations from shared effect by |posterior mean| ===\n")
print(head(top_deviations, 20), row.names = FALSE)

# -- 8. Write CSVs -------------------------------------------------------------
out_files <- c(
  scan_full      = paste0(OUT_PREFIX, "_ci_scan_full_vs.csv"),
  beta80         = paste0(OUT_PREFIX, "_beta_80_included_vs.csv"),
  beta90         = paste0(OUT_PREFIX, "_beta_90_included_vs.csv"),
  beta95         = paste0(OUT_PREFIX, "_beta_95_included_vs.csv"),
  delta80        = paste0(OUT_PREFIX, "_delta_80_included_vs.csv"),
  delta90        = paste0(OUT_PREFIX, "_delta_90_included_vs.csv"),
  delta95        = paste0(OUT_PREFIX, "_delta_95_included_vs.csv"),
  gene_summary   = paste0(OUT_PREFIX, "_gene_summary_vs.csv"),
  ezh2_kmt2d     = paste0(OUT_PREFIX, "_ezh2_kmt2d_vs.csv"),
  braf_egfr      = paste0(OUT_PREFIX, "_braf_egfr_vs.csv"),
  top_beta       = paste0(OUT_PREFIX, "_top_beta_effects_vs.csv"),
  top_deviation  = paste0(OUT_PREFIX, "_top_deviations_vs.csv")
)

write.csv(scan_df, out_files[["scan_full"]], row.names = FALSE)
write.csv(beta80, out_files[["beta80"]], row.names = FALSE)
write.csv(beta90, out_files[["beta90"]], row.names = FALSE)
write.csv(beta95, out_files[["beta95"]], row.names = FALSE)
write.csv(delta80, out_files[["delta80"]], row.names = FALSE)
write.csv(delta90, out_files[["delta90"]], row.names = FALSE)
write.csv(delta95, out_files[["delta95"]], row.names = FALSE)
write.csv(gene_summary, out_files[["gene_summary"]], row.names = FALSE)
write.csv(ezh2_kmt2d, out_files[["ezh2_kmt2d"]], row.names = FALSE)
write.csv(braf_egfr, out_files[["braf_egfr"]], row.names = FALSE)
write.csv(top_beta_effects, out_files[["top_beta"]], row.names = FALSE)
write.csv(top_deviations, out_files[["top_deviation"]], row.names = FALSE)

cat("\n=== FILES WRITTEN ===\n")
for (f in out_files) {
  if (file.exists(f)) cat(sprintf("  [OK] %s (%.1f KB)\n", f, file.size(f) / 1024))
}

cat("\n=== SUMMARY FOR MANUSCRIPT / NOTES ===\n")
cat(sprintf("Total tissue-target-gene-drug combinations scanned: %d\n", nrow(scan_df)))
cat(sprintf("Target genes scanned: %d\n", T_TARGET))
cat(sprintf("Groups scanned: %d\n", K))
cat(sprintf("Drugs scanned: %d\n", D))
cat(sprintf("Beta 80%% CI flags: %d\n", nrow(beta80)))
cat(sprintf("Beta 90%% CI flags: %d\n", nrow(beta90)))
cat(sprintf("Beta 95%% CI flags: %d\n", nrow(beta95)))
cat(sprintf("Delta 80%% CI flags: %d\n", nrow(delta80)))
cat(sprintf("Delta 90%% CI flags: %d\n", nrow(delta90)))
cat(sprintf("Delta 95%% CI flags: %d\n", nrow(delta95)))
cat(sprintf("EZH2/KMT2D rows scanned: %d\n", nrow(ezh2_kmt2d)))
cat(sprintf("BRAF/EGFR/TP53/ASXL1/NRAS/STAG2 rows scanned: %d\n", nrow(braf_egfr)))
cat("\nDone.\n")
