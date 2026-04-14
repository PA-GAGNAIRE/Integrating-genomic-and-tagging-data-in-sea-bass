################################################################################
# Figure 3 – Individual ancestry proportions from ADMIXTURE
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 3 ##
#
# CHANGES FROM RMD (for standalone reusability):
#   [1] Added library() calls at the top of the script.
#   [2] Added data-loading block from the "Data loading" chunk in the Rmd.
#   [3] Added configurable path variables (data_dir, figures_dir) so users
#       can adapt the script to their own directory layout without editing
#       the core code.
#   All plotting code below is UNCHANGED from the Rmd.
################################################################################

# [CHANGE 1] Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
metadata  <- read.csv(file.path(data_dir, "metadata.csv"))
admix_K2  <- read.table(file.path(data_dir, "admixture_K2.txt"))
admix_K3  <- read.table(file.path(data_dir, "admixture_K3.txt"))

# Create figures directory if it does not already exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- UNCHANGED CODE FROM RMD BELOW ----------

# Attach metadata to admixture results (individuals must be in the same order)
admix_K2_df <- cbind(metadata, admix_K2)
admix_K3_df <- cbind(metadata, admix_K3)

# Geographic order: sort populations from west (Atlantic) to east (Mediterranean)
geo_order <- metadata %>%
  group_by(site, region) %>%
  summarise(mean_lon = mean(longitude), .groups = "drop") %>%
  arrange(mean_lon) %>%
  pull(site)

# Helper function to reshape and plot admixture bar charts
plot_admixture <- function(df, K, colors, title) {
  q_cols <- paste0("V", seq_len(K))
  df_long <- df %>%
    mutate(individual_id = factor(individual_id,
                                  levels = df$individual_id[order(match(df$site, geo_order),
                                                                  df[[q_cols[1]]])]),
           site = factor(site, levels = geo_order)) %>%
    pivot_longer(cols = all_of(q_cols),
                 names_to  = "cluster",
                 values_to = "ancestry")

  ggplot(df_long, aes(x = individual_id, y = ancestry, fill = cluster)) +
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(values = colors, name = "Cluster") +
    facet_grid(~ site, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Ancestry proportion", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.text       = element_text(angle = 45, hjust = 1, size = 8),
      panel.spacing    = unit(0.05, "lines"),
      legend.position  = "bottom"
    )
}

p_K2 <- plot_admixture(admix_K2_df, K = 2,
                       colors = c("#2166AC", "#D73027"), title = "K = 2")
p_K3 <- plot_admixture(admix_K3_df, K = 3,
                       colors = c("#2166AC", "#4DAC26", "#D73027"), title = "K = 3")

figure3 <- p_K2 / p_K3 +
  plot_annotation(tag_levels = "A")

figure3

ggsave(file.path(figures_dir, "Figure3_admixture.pdf"), figure3,
       width = 14, height = 6, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure3_admixture.png"), figure3,
       width = 14, height = 6, dpi = 300)
