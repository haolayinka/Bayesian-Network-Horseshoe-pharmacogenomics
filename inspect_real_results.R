#!/usr/bin/env Rscript
#
# inspect_real_results.R
#
# Inspects the posterior results from a completed full production run
# (fits_full), produced by gibbs_sampler_full.Rmd Part 3 Step 4. This is
# the first real scientific look at the model's output: which genes show
# up as drivers of drug response, starting with TP53 (the panel's most
# frequently mutated gene, with an external precedent in Samorodnitsky et
# al. 2020 for survival), then a full gene x drug inclusion scan, then a
# specific check of the MSI co-mutation cluster genes already flagged in
# the multicollinearity diagnostic.
#
# CRITICAL DEPENDENCY -- READ BEFORE RUNNING:
# gibbs_sampler()'s output (fit$beta_tilde, an array of dimension
# [G, D, n_kept]) does NOT store which drug_id corresponds to which
# column index in the D dimension -- this information exists only in
# however y_real's columns were ordered when fits_full was created. In
# gibbs_sampler_full.Rmd, y_real's columns are built from:
#     all_drugs <- sort(unique(resp$drug_id))
# in that exact sorted order. The reconstruction below replicates this
# EXACTLY. If you used a different drug ordering or a subset when
# generating fits_full, you MUST adjust `all_drugs_used` below to match,
# or every drug-level result in this script will be silently mislabeled
# -- this was verified to be a real risk during development (an earlier
# test run on a 30-drug subset confirmed the labeling is correct only
# when this reconstruction matches the original ordering exactly).
#
# Usage: edit Step 1 below to point at your saved fits_full.rds and your
# real data files, then run top to bottom (or source this from an
# interactive R session after loading fits_full yourself).

# --------------------------------------------------------------------- #
# Step 1: Load your completed fit and real data
# --------------------------------------------------------------------- #

fits_full <- readRDS("gibbs_fits_full.rds")
mut <- read.csv("mutation_matrix.csv", row.names = 1, check.names = FALSE)
resp <- read.csv("response_matrix.csv")

gene_names <- fits_full[[1]]$gene_names
n_genes <- length(gene_names)
n_drugs <- dim(fits_full[[1]]$beta_tilde)[2]
n_chains <- length(fits_full)

# Reconstruct drug IDs in the SAME order used to build y_real in
# gibbs_sampler_full.Rmd (sort(unique(resp$drug_id))). If fits_full was
# built on the FULL drug panel (the normal case for a production run),
# this is all 295 drugs and no further editing is needed below. If you
# fit a SUBSET of drugs, edit all_drugs_used to match that subset exactly
# (e.g. all_drugs_used <- sort(unique(resp$drug_id))[1:30]).
all_drugs_used <- sort(unique(resp$drug_id))

stopifnot(length(all_drugs_used) == n_drugs)  # hard stop if the reconstruction doesn't match
cat(sprintf("Loaded %d chains, %d genes, %d drugs.\n", n_chains, n_genes, n_drugs))
cat("First 5 reconstructed drug IDs:", head(all_drugs_used, 5), "\n")
cat("Last 5 reconstructed drug IDs:", tail(all_drugs_used, 5), "\n")

# --------------------------------------------------------------------- #
# Step 2: Helper to pool posterior samples across all chains for one
# gene x drug pair, and compute a 95% credible interval
# --------------------------------------------------------------------- #

get_pooled_ci <- function(gene_idx, drug_idx, prob = c(0.025, 0.975)) {
  vals <- unlist(lapply(fits_full, function(f) f$beta_tilde[gene_idx, drug_idx, ]))
  ci <- quantile(vals, prob)
  list(post_mean = mean(vals), post_sd = sd(vals), ci_lo = ci[1], ci_hi = ci[2],
       excludes_zero = (ci[1] > 0) | (ci[2] < 0))
}

# --------------------------------------------------------------------- #
# Step 3: TP53 across all drugs
# --------------------------------------------------------------------- #

tp53_idx <- which(gene_names == "TP53")
if (length(tp53_idx) == 0) stop("TP53 not found in gene_names -- check the gene panel.")

tp53_summary <- data.frame(drug_id = all_drugs_used, post_mean = NA_real_,
                            post_sd = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                            excludes_zero = NA)
for (d in seq_len(n_drugs)) {
  r <- get_pooled_ci(tp53_idx, d)
  tp53_summary$post_mean[d] <- r$post_mean
  tp53_summary$post_sd[d] <- r$post_sd
  tp53_summary$ci_lo[d] <- r$ci_lo
  tp53_summary$ci_hi[d] <- r$ci_hi
  tp53_summary$excludes_zero[d] <- r$excludes_zero
}

cat(sprintf("\n=== TP53 across all %d drugs ===\n", n_drugs))
cat("Drugs where TP53's 95%% CI excludes zero:", sum(tp53_summary$excludes_zero), "/", n_drugs, "\n")
cat("\nTop 10 drugs by |posterior mean| (strongest TP53 effect, either direction):\n")
print(head(tp53_summary[order(-abs(tp53_summary$post_mean)), ], 10))

write.csv(tp53_summary, "tp53_across_drugs.csv", row.names = FALSE)

# --------------------------------------------------------------------- #
# Step 4: Full inclusion scan -- every gene x every drug
#
# WARNING: this is G x D credible-interval computations (e.g. 219 x 295
# = 64,605 for the full real panel) -- expect this to take a few minutes,
# not seconds, at full scale. Consider running Step 3 (TP53 only) and
# Step 5 (MSI cluster only) first if you want fast initial results before
# committing to the full scan.
# --------------------------------------------------------------------- #

cat(sprintf("\n=== Full inclusion scan: %d genes x %d drugs = %d pairs ===\n",
            n_genes, n_drugs, n_genes * n_drugs))
cat("This may take several minutes at full scale...\n")

inclusion_scan <- vector("list", n_genes * n_drugs)
idx <- 1L
for (g in seq_len(n_genes)) {
  for (d in seq_len(n_drugs)) {
    r <- get_pooled_ci(g, d)
    inclusion_scan[[idx]] <- data.frame(
      gene = gene_names[g], drug_id = all_drugs_used[d],
      post_mean = r$post_mean, post_sd = r$post_sd,
      ci_lo = r$ci_lo, ci_hi = r$ci_hi, excludes_zero = r$excludes_zero
    )
    idx <- idx + 1L
  }
}
inclusion_scan <- do.call(rbind, inclusion_scan)

cat("Pairs where 95% CI excludes zero:", sum(inclusion_scan$excludes_zero),
    "/", nrow(inclusion_scan), "\n\n")

included <- inclusion_scan[inclusion_scan$excludes_zero, ]
gene_counts <- sort(table(included$gene), decreasing = TRUE)
cat("Genes appearing as a 'driver' (CI excludes zero) for at least one drug,\n")
cat("ranked by number of drugs (top 20):\n")
print(head(gene_counts, 20))

write.csv(inclusion_scan, "full_inclusion_scan.csv", row.names = FALSE)
write.csv(included, "included_gene_drug_pairs.csv", row.names = FALSE)
cat("\nSaved full_inclusion_scan.csv and included_gene_drug_pairs.csv\n")

# --------------------------------------------------------------------- #
# Step 5: MSI co-mutation cluster check
#
# RPL22, ACVR2A, BAX, BMPR2, KMT2B, STAT5B were flagged in the
# multicollinearity diagnostic (condition number 208.12) as a tight,
# biologically-explained co-mutation cluster linked to microsatellite
# instability (MSI). Worth checking directly whether this model finds
# coordinated "driver" signal across the same drugs for these genes,
# which would itself be a reportable finding tied to that diagnostic.
# --------------------------------------------------------------------- #

msi_cluster <- c("RPL22", "ACVR2A", "BAX", "BMPR2", "KMT2B", "STAT5B")
msi_present <- intersect(msi_cluster, gene_names)
cat(sprintf("\n=== MSI cluster check (%d of %d cluster genes present in panel) ===\n",
            length(msi_present), length(msi_cluster)))

msi_scan <- inclusion_scan[inclusion_scan$gene %in% msi_present, ]
cat("MSI-cluster gene-drug pairs where CI excludes zero:",
    sum(msi_scan$excludes_zero), "/", nrow(msi_scan), "\n\n")

msi_included <- msi_scan[msi_scan$excludes_zero, ]
if (nrow(msi_included) > 0) {
  cat("Drugs where an MSI-cluster gene shows a driver signal:\n")
  print(table(msi_included$drug_id, msi_included$gene))

  # Check for coordinated signal: do multiple MSI genes show up as drivers
  # for the SAME drug? This is the specific pattern worth flagging if
  # present, since it would suggest the model is picking up the shared MSI
  # mechanism rather than (or in addition to) gene-specific effects.
  drugs_with_multiple_msi_hits <- table(msi_included$drug_id)
  drugs_with_multiple_msi_hits <- drugs_with_multiple_msi_hits[drugs_with_multiple_msi_hits > 1]
  if (length(drugs_with_multiple_msi_hits) > 0) {
    cat("\nDrugs where >=2 MSI-cluster genes BOTH show a driver signal\n")
    cat("(possible coordinated MSI-mechanism signal, not just gene-specific):\n")
    print(drugs_with_multiple_msi_hits)
  } else {
    cat("\nNo drug shows >=2 MSI-cluster genes simultaneously as drivers --\n")
    cat("any signal found is gene-specific, not obviously a coordinated MSI effect.\n")
  }
} else {
  cat("No MSI-cluster gene shows a driver signal (CI excludes zero) for any drug.\n")
}

write.csv(msi_scan, "msi_cluster_scan.csv", row.names = FALSE)
cat("\nSaved msi_cluster_scan.csv\n")
