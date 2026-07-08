# ============================================================
# fig2_rmse_decomposed.R
#
# Generates Figure: RMSE decomposed by true-association status
# (true-nonzero signal pairs vs true-zero null pairs)
#
# UPDATE SCENARIO 2 VALUES HERE when new simulation run completes.
# All other values are final (n=100 replicates).
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Data: update sc2 values when Scenario 2 rerun completes ─────────────
rmse_data <- tribble(
  ~scenario,              ~model,         ~signal, ~null,
  # Scenario 1 (FINAL - n=100)
  "Sc.1: Sparse\nIndependent", "Full model",    0.562,  0.066,
  "Sc.1: Sparse\nIndependent", "No-network",    0.263,  0.102,
  "Sc.1: Sparse\nIndependent", "No-horseshoe",  0.415,  0.166,
  # Scenario 2 (UPDATE THESE when new run completes)
  "Sc.2: Network\nStructured", "Full model",    0.586,  0.067,
  "Sc.2: Network\nStructured", "No-network",    0.263,  0.099,
  "Sc.2: Network\nStructured", "No-horseshoe",  0.418,  0.166,
  # Scenario 3 (FINAL - n=100)
  "Sc.3: Dense\nWeak Noise",   "Full model",    0.562,  0.073,
  "Sc.3: Dense\nWeak Noise",   "No-network",    0.259,  0.104,
  "Sc.3: Dense\nWeak Noise",   "No-horseshoe",  0.411,  0.167
)

# ── Reshape to long format ───────────────────────────────────────────────
rmse_long <- rmse_data %>%
  pivot_longer(c(signal, null),
               names_to  = "subset",
               values_to = "rmse") %>%
  mutate(
    scenario = factor(scenario,
                      levels = c("Sc.1: Sparse\nIndependent",
                                 "Sc.2: Network\nStructured",
                                 "Sc.3: Dense\nWeak Noise")),
    model    = factor(model,
                      levels = c("Full model", "No-network", "No-horseshoe")),
    subset   = factor(subset,
                      levels = c("signal", "null"),
                      labels = c("True nonzero pairs\n(3.3% of all pairs)",
                                 "True zero pairs\n(96.7% of all pairs)"))
  )

# ── Colours (colour-blind friendly) ─────────────────────────────────────
model_colours <- c(
  "Full model"    = "#0072B2",
  "No-network"    = "#E69F00",
  "No-horseshoe"  = "#CC79A7"
)

# ── Plot ─────────────────────────────────────────────────────────────────
p <- ggplot(rmse_long,
            aes(x      = model,
                y      = rmse,
                fill   = model,
                alpha  = subset,
                pattern = subset)) +

  # Bars -- solid for signal, hatched for null
  geom_col(position = position_dodge(width = 0.75),
           width    = 0.35,
           colour   = "white", linewidth = 0.3) +

  # Value labels
  geom_text(aes(label    = sprintf("%.3f", rmse),
                fontface = ifelse(subset ==
                                    "True nonzero pairs\n(3.3% of all pairs)",
                                  "bold", "plain")),
            position = position_dodge(width = 0.75),
            vjust    = -0.4,
            size     = 2.8) +

  facet_wrap(~scenario, nrow = 1) +

  scale_fill_manual(values = model_colours, name = NULL) +
  scale_alpha_manual(
    values = c("True nonzero pairs\n(3.3% of all pairs)" = 0.87,
               "True zero pairs\n(96.7% of all pairs)"   = 0.40),
    name   = NULL
  ) +

  scale_y_continuous(limits = c(0, 0.72),
                     breaks = seq(0, 0.7, by = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +

  labs(
    title = "Figure 2  |  RMSE Decomposed by True-Association Status (n = 100 replicates)",
    x     = NULL,
    y     = "RMSE"
  ) +

  theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.key.size   = unit(0.5, "cm"),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background  = element_rect(fill = "grey92"),
    strip.text        = element_text(face = "bold", size = 10),
    plot.title        = element_text(face = "bold", size = 11),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank()
  ) +

  # Separate legends for model colour and subset pattern
  guides(
    fill  = guide_legend(order = 2, nrow = 1,
                         override.aes = list(alpha = 0.87)),
    alpha = guide_legend(order = 1, nrow = 1,
                         override.aes = list(fill = "grey50"))
  )

# ── Save at 300 DPI ──────────────────────────────────────────────────────
ggsave("fig2_rmse_decomposed.png",
       plot   = p,
       width  = 13,
       height = 5,
       dpi    = 300,
       units  = "in",
       bg     = "white")

cat("Saved: fig2_rmse_decomposed.png\n")
cat("\nTo update Scenario 2: edit the three sc2 rows in rmse_data above\n")
cat("and re-run the script. Everything else is automatic.\n")
