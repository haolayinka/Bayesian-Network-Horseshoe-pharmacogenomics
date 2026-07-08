# R Code: Network-Structured Bayesian Hierarchical Model for Sparse Mutation–Drug Response Associations

Code accompanying the manuscript:
**"A Network-Structured Bayesian Hierarchical Model for Sparse Mutation–Drug Response Associations: Application to Cancer Pharmacogenomics"**

---

## File Descriptions

### Core Model

| File | Description |
|------|-------------|
| `gibbs_sampler.R` | Primary Gibbs sampler for Model I (shared-effect network-horseshoe model). Implements the GIGG conjugate updates for all parameters. |
| `gibbs_sampler_full_AB.Rmd` | RMarkdown version of the primary sampler with full inline documentation and derivations. |
| `real_data_group_layer_vs_target_storage.R` | Gibbs sampler for Model II (tissue-group extension). Implements the group-specific coefficient updates with IG(2,0.5) prior on eta^2 and ridge stabilization. Minimum group size = 50. |

### Real Data Analysis

| File | Description |
|------|-------------|
| `inspect_real_results.R` | Extracts and summarises primary GDSC2 results from Model I fit objects. |
| `scan_tissue_effects_vs_target_posterior.R` | Extracts and summarises tissue-group results from Model II fit objects (K=6). |
| `inspect_prism_results.R` | Extracts and summarises PRISM/DepMap replication results. |
| `prism_pipeline.R` | Preprocesses PRISM data for compatibility with the model framework. |

### Validation and Sensitivity

| File | Description |
|------|-------------|
| `kfold_cv_comparison.Rmd` | Five-fold cross-validated predictive log-likelihood comparison across all three model configurations (full, no-network, no-horseshoe). |
| `network_misspecification_check.R` | Delta sensitivity analysis: reruns the model at delta in {0.5, 0.7, 0.8, 0.9} and compares flagged associations. |
| `run_simulation_test.R` | Runs simulation study replicates under three data-generating scenarios. |

### Figures and Tables

| File | Description |
|------|-------------|
| `fig2_rmse_decomposed.R` | Generates Figure 2: RMSE decomposed into signal and null subsets across scenarios. |
| `fig2_rmse_decomposed.py` | Python equivalent of Figure 2 script. |
| `rmse_decomposition_table.Rmd` | Generates the RMSE decomposition table (requires pair_scores/ folder). |
| `simulation_plots.Rmd` | Generates simulation study figures. |

---

## Data Requirements

The following data files are required and should be placed in the working directory:

- `mutation_matrix.csv` — binary somatic mutation matrix (samples × genes)
- `response_matrix.csv` — drug response matrix with columns: model_id, drug_id, drug_name, y
- `model_groups.csv` — tissue-of-origin annotations with columns: model_id, tissue
- `gene_adjacency_W.csv` — 135-gene pathway adjacency matrix (derived from KEGG)
- `pathway_gene_membership.csv` — KEGG pathway gene membership (used to construct adjacency)

GDSC2 data: https://www.cancerrxgene.org  
PRISM data: https://depmap.org/portal  
KEGG pathway data: https://www.genome.jp/kegg

---

## Key Implementation Notes

1. **GIGG:** Both global variance components (kappa^2, zeta^2) are parameterized through their precisions with Gamma(1,1) priors. Do NOT revert to IG(0.01,0.01) on the variances directly — this causes both to collapse to zero.

2. **Adjacency matrix:** Always use `gene_adjacency_W.csv` (gene-level, 135 nodes). Do not use pathway-level adjacency.

3. **Group layer minimum size:** Set MIN_GROUP_SIZE = 50 in `real_data_group_layer_v2.R`. Groups below this threshold are merged into "Other".

4. **MCMC settings:** 4,000 iterations, 1,000 burn-in, thinning = 2, 2 chains per model configuration.

---

## Session Info

R version 4.x; key packages: `Matrix`, `dplyr`, `ggplot2`, `coda`
