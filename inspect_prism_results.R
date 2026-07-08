#!/usr/bin/env Rscript
#
# inspect_prism_results.R
#
# PRISM/DepMap external-validation version of inspect_real_results.R.
# Inspects the posterior results from the PRISM Repurposing 24Q2 run
# (fits produced by running gibbs_sampler.R on the PRISM pipeline outputs).
#
# KEY DIFFERENCES FROM inspect_real_results.R:
#   1. Points at PRISM files (mutation_matrix_prism.csv, response_matrix_prism.csv)
#   2. drug_id is a broad_id string (e.g. BRD-K12345-001-01-1), not numeric
#   3. Adds drug_name lookup from response matrix so output is human-readable
#   4. Primary check is EZH2 and KMT2D (not TP53) -- the GDSC2 headline findings
#   5. G=218 genes (one fewer than GDSC2), D=1518 compounds
#
# CRITICAL DEPENDENCY (same as GDSC2 version):
#   all_drugs_used MUST be reconstructed as sort(unique(resp$drug_id)) to
#   match the column ordering of y_real used when fits were created.
#   The stopifnot() below will catch any mismatch immediately.
#
# Usage: edit Step 1 to point at your saved .rds files, then run top to bottom.

# --------------------------------------------------------------------- #
# Step 1: Load completed fits and PRISM data files
# --------------------------------------------------------------------- #

fits_full       <- readRDS("gibbs_fits_full_prism.rds")       # adjust path as needed
fits_nonetwork  <- readRDS("gibbs_fits_nonetwork_prism.rds")  # adjust path as needed
fits_nohorseshoe <- readRDS("gibbs_fits_nohorseshoe_prism.rds") # adjust path as needed

mut  <- read.csv("mutation_matrix_prism.csv",  row.names = 1, check.names = FALSE)
resp <- read.csv("response_matrix_prism.csv")

cat("Response matrix columns:", paste(colnames(resp), collapse=", "), "\n")

# Build drug name lookup from the response matrix
# (drug_name column was added by prism_pipeline.R from PortalCompounds.csv)
drug_lookup <- unique(resp[, c("drug_id", "drug_name")])
drug_lookup <- drug_lookup[order(drug_lookup$drug_id), ]
rownames(drug_lookup) <- NULL

get_drug_name <- function(drug_id) {
  idx <- match(drug_id, drug_lookup$drug_id)
  ifelse(is.na(idx), drug_id, drug_lookup$drug_name[idx])
}

gene_names <- fits_full[[1]]$gene_names
n_genes    <- length(gene_names)
n_drugs    <- dim(fits_full[[1]]$beta_tilde)[2]
n_chains   <- length(fits_full)

# Reconstruct drug IDs in EXACT same order as y_real was built
all_drugs_used <- sort(unique(resp$drug_id))

# Hard stop if reconstruction doesn't match -- catches any mismatch immediately
stopifnot(length(all_drugs_used) == n_drugs)

cat(sprintf("Loaded: %d chains | %d genes | %d compounds\n",
            n_chains, n_genes, n_drugs))
cat("First 3 drug IDs:  ", paste(head(all_drugs_used, 3), collapse=", "), "\n")
cat("First 3 drug names:", paste(get_drug_name(head(all_drugs_used, 3)), collapse=", "), "\n")

# --------------------------------------------------------------------- #
# Step 2: Helper -- pool posterior samples across all chains for one
#         gene x drug pair
# --------------------------------------------------------------------- #

get_pooled_ci <- function(fits, gene_idx, drug_idx, prob = c(0.025, 0.975)) {
  vals <- unlist(lapply(fits, function(f) f$beta_tilde[gene_idx, drug_idx, ]))
  ci   <- quantile(vals, prob)
  list(post_mean     = mean(vals),
       post_sd       = sd(vals),
       ci_lo         = ci[1],
       ci_hi         = ci[2],
       excludes_zero = (ci[1] > 0) | (ci[2] < 0))
}

# --------------------------------------------------------------------- #
# Step 3: PRIMARY CHECK -- EZH2 and KMT2D across all compounds
#         These are the GDSC2 headline findings; replication here is the
#         main scientific contribution of the PRISM validation.
# --------------------------------------------------------------------- #

cat("\n=== PRIMARY CHECK: EZH2 and KMT2D across all compounds ===\n")

for (target_gene in c("EZH2", "KMT2D")) {
  g_idx <- which(gene_names == target_gene)
  if (length(g_idx) == 0) {
    cat(sprintf("WARNING: %s not found in gene panel -- skipped\n", target_gene))
    next
  }

  gene_summary <- data.frame(
    drug_id       = all_drugs_used,
    drug_name     = get_drug_name(all_drugs_used),
    post_mean     = NA_real_,
    post_sd       = NA_real_,
    ci_lo         = NA_real_,
    ci_hi         = NA_real_,
    excludes_zero = NA
  )

  for (d in seq_len(n_drugs)) {
    r <- get_pooled_ci(fits_full, g_idx, d)
    gene_summary$post_mean[d]     <- r$post_mean
    gene_summary$post_sd[d]       <- r$post_sd
    gene_summary$ci_lo[d]         <- r$ci_lo
    gene_summary$ci_hi[d]         <- r$ci_hi
    gene_summary$excludes_zero[d] <- r$excludes_zero
  }

  n_flagged   <- sum(gene_summary$excludes_zero)
  n_negative  <- sum(gene_summary$post_mean[gene_summary$excludes_zero] < 0)
  mean_effect <- mean(gene_summary$post_mean[gene_summary$excludes_zero])

  cat(sprintf("\n%s: %d / %d compounds flagged (95%% CI excludes zero)\n",
              target_gene, n_flagged, n_drugs))
  if (n_flagged > 0) {
    cat(sprintf("  Direction: %d negative (sensitivity) | %d positive (resistance)\n",
                n_negative, n_flagged - n_negative))
    cat(sprintf("  Mean posterior effect: %.4f (GDSC2 was %.3f)\n",
                mean_effect,
                ifelse(target_gene == "EZH2", -0.911, -0.496)))
    cat("  Top 10 by |posterior mean|:\n")
    sub <- gene_summary[gene_summary$excludes_zero, ]
    sub <- sub[order(-abs(sub$post_mean)), ]
    print(head(sub[, c("drug_name","post_mean","post_sd","ci_lo","ci_hi")], 10),
          row.names = FALSE)
  }

  fname <- sprintf("%s_prism_scan.csv", tolower(target_gene))
  write.csv(gene_summary, fname, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", fname))
}

# --------------------------------------------------------------------- #
# Step 4: TP53 -- check for Nutlin-3a replication (positive control)
# --------------------------------------------------------------------- #

cat("\n=== TP53 check (Nutlin-3a positive control) ===\n")
tp53_idx <- which(gene_names == "TP53")

if (length(tp53_idx) > 0) {
  # Check if Nutlin-3a is in the PRISM compound panel
  nutlin_idx <- grep("nutlin", drug_lookup$drug_name, ignore.case = TRUE)
  if (length(nutlin_idx) > 0) {
    cat("Nutlin-3a found in PRISM panel:\n")
    print(drug_lookup[nutlin_idx, ], row.names = FALSE)
    # Check TP53 effect on Nutlin-3a
    for (ni in nutlin_idx) {
      d_idx <- which(all_drugs_used == drug_lookup$drug_id[ni])
      if (length(d_idx) > 0) {
        r <- get_pooled_ci(fits_full, tp53_idx, d_idx)
        cat(sprintf("  TP53 effect on %s: mean=%.4f, SD=%.4f, CI=[%.4f, %.4f], excl_zero=%s\n",
                    drug_lookup$drug_name[ni],
                    r$post_mean, r$post_sd, r$ci_lo, r$ci_hi,
                    r$excludes_zero))
        cat(sprintf("  GDSC2 reference: mean=+1.351, CI=[1.063, 1.612] (resistance direction)\n"))
      }
    }
  } else {
    cat("Nutlin-3a NOT found in PRISM 1,518-compound panel.\n")
    cat("(PRISM Repurposing focuses on approved/clinical-stage compounds;\n")
    cat(" Nutlin-3a is a tool compound and may not be included.)\n")
  }
} else {
  cat("TP53 not in gene panel.\n")
}

# --------------------------------------------------------------------- #
# Step 5: Full inclusion scan -- all 218 x 1518 pairs
#
# WARNING: 218 x 1518 = 330,924 pairs -- expect 15-30 minutes at full scale.
# Steps 3-4 above give the headline results quickly.
# Comment this section out if you only need the primary check.
# --------------------------------------------------------------------- #

cat(sprintf("\n=== Full inclusion scan: %d genes x %d compounds = %d pairs ===\n",
            n_genes, n_drugs, n_genes * n_drugs))
cat("This may take 15-30 minutes. Steps 3-4 above have the headline results.\n")
cat("Running...\n")

inclusion_scan <- vector("list", n_genes * n_drugs)
idx <- 1L
for (g in seq_len(n_genes)) {
  for (d in seq_len(n_drugs)) {
    r <- get_pooled_ci(fits_full, g, d)
    inclusion_scan[[idx]] <- data.frame(
      gene          = gene_names[g],
      drug_id       = all_drugs_used[d],
      drug_name     = get_drug_name(all_drugs_used[d]),
      post_mean     = r$post_mean,
      post_sd       = r$post_sd,
      ci_lo         = r$ci_lo,
      ci_hi         = r$ci_hi,
      excludes_zero = r$excludes_zero
    )
    idx <- idx + 1L
  }
}
inclusion_scan <- do.call(rbind, inclusion_scan)

cat(sprintf("Pairs where 95%% CI excludes zero: %d / %d (%.3f%%)\n",
            sum(inclusion_scan$excludes_zero),
            nrow(inclusion_scan),
            100 * mean(inclusion_scan$excludes_zero)))

included <- inclusion_scan[inclusion_scan$excludes_zero, ]
cat("\nTop driver genes (ranked by number of associated compounds):\n")
print(head(sort(table(included$gene), decreasing = TRUE), 20))

cat("\n=== EZH2 and KMT2D direction check (all flagged compounds) ===\n")
for (g in c("EZH2", "KMT2D")) {
  sub <- included[included$gene == g, ]
  if (nrow(sub) > 0) {
    cat(sprintf("%s: %d compounds | all negative: %s | mean effect: %.4f\n",
                g, nrow(sub),
                all(sub$post_mean < 0),
                mean(sub$post_mean)))
  } else {
    cat(sprintf("%s: 0 compounds flagged\n", g))
  }
}

write.csv(inclusion_scan, "full_inclusion_scan_prism.csv",   row.names = FALSE)
write.csv(included,       "included_gene_drug_pairs_prism.csv", row.names = FALSE)
cat("\nSaved: full_inclusion_scan_prism.csv and included_gene_drug_pairs_prism.csv\n")

# --------------------------------------------------------------------- #
# Step 6: Ablation comparison
# --------------------------------------------------------------------- #

cat("\n=== Ablation comparison ===\n")

count_flagged <- function(fits) {
  total <- 0L
  for (g in seq_len(n_genes))
    for (d in seq_len(n_drugs)) {
      r <- get_pooled_ci(fits, g, d)
      if (r$excludes_zero) total <- total + 1L
    }
  total
}

cat("Counting flagged pairs for no-network model...\n")
n_nonetwork    <- count_flagged(fits_nonetwork)
cat("Counting flagged pairs for no-horseshoe model...\n")
n_nohorseshoe  <- count_flagged(fits_nohorseshoe)
n_full         <- sum(inclusion_scan$excludes_zero)

cat(sprintf("\nFull model:     %d flagged (%.3f%%)\n", n_full,
            100*n_full/(n_genes*n_drugs)))
cat(sprintf("No-network:     %d flagged (%.3f%%) | ratio to full: %.1fx\n",
            n_nonetwork, 100*n_nonetwork/(n_genes*n_drugs),
            n_nonetwork/max(n_full,1)))
cat(sprintf("No-horseshoe:   %d flagged (%.3f%%) | ratio to full: %.1fx\n",
            n_nohorseshoe, 100*n_nohorseshoe/(n_genes*n_drugs),
            n_nohorseshoe/max(n_full,1)))
cat(sprintf("Ordering (full < no-network < no-horseshoe): %s\n",
            n_full < n_nonetwork && n_nonetwork < n_nohorseshoe))

cat("\n=== PRISM inspection complete ===\n")
