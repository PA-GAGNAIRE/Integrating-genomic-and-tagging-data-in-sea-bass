################################################################################
# Figure 4 – Pairwise FST heatmap and Isolation-by-Distance (IBD) test
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 4 ##
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
library(tibble)
library(patchwork)
library(vegan)
library(fields)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
metadata   <- read.csv(file.path(data_dir, "metadata.csv"))
fst_matrix <- read.csv(file.path(data_dir, "fst_matrix.csv"), row.names = 1)

# Create figures directory if it does not already exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- UNCHANGED CODE FROM RMD BELOW ----------

# Geographic order: sort populations from west (Atlantic) to east (Mediterranean)
geo_order <- metadata %>%
  group_by(site, region) %>%
  summarise(mean_lon = mean(longitude), .groups = "drop") %>%
  arrange(mean_lon) %>%
  pull(site)

# --- FST heatmap ---
fst_long <- fst_matrix %>%
  rownames_to_column("Pop1") %>%
  pivot_longer(-Pop1, names_to = "Pop2", values_to = "Fst") %>%
  mutate(Fst = ifelse(is.na(Fst) | Pop1 == Pop2, NA_real_, Fst))

# Order populations geographically (west to east)
pop_order <- geo_order[geo_order %in% rownames(fst_matrix)]
fst_long$Pop1 <- factor(fst_long$Pop1, levels = pop_order)
fst_long$Pop2 <- factor(fst_long$Pop2, levels = pop_order)

p_fst_heat <- ggplot(fst_long, aes(x = Pop2, y = Pop1, fill = Fst)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(!is.na(Fst), sprintf("%.3f", Fst), "")),
            size = 2.8) +
  scale_fill_gradient2(low = "white", mid = "#FEE090", high = "#D73027",
                       midpoint = median(fst_long$Fst, na.rm = TRUE),
                       na.value = "grey95", name = expression(italic(F)[ST])) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# --- Isolation by Distance ---
# Compute mean coordinates per population
pop_coords <- metadata %>%
  group_by(site) %>%
  summarise(lon = mean(longitude), lat = mean(latitude), .groups = "drop") %>%
  filter(site %in% pop_order)

# Geographic distance matrix (great-circle, in km)
coord_mat    <- as.matrix(pop_coords[, c("lon", "lat")])
rownames(coord_mat) <- pop_coords$site
geo_dist     <- as.dist(fields::rdist.earth(coord_mat, miles = FALSE))

# FST / (1 - FST) linearization
fst_sym      <- as.matrix(fst_matrix)[pop_order, pop_order]
fst_lin      <- as.dist(fst_sym / (1 - fst_sym))

# Mantel test
mantel_res   <- vegan::mantel(fst_lin, geo_dist, method = "pearson",
                               permutations = 999)

# Data frame for scatter plot
ibd_df <- data.frame(
  geo_dist = as.vector(geo_dist) / 1000,   # km -> 1000 km
  fst_lin  = as.vector(fst_lin)
)

p_ibd <- ggplot(ibd_df, aes(x = geo_dist, y = fst_lin)) +
  geom_point(alpha = 0.6, color = "#2166AC", size = 2.5) +
  geom_smooth(method = "lm", color = "#D73027", fill = "#FDAE61",
              linewidth = 1, se = TRUE) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Mantel r = %.3f\np = %.3f",
                           mantel_res$statistic, mantel_res$signif),
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(x = "Geographic distance (× 1,000 km)",
       y = expression(italic(F)[ST] / (1 - italic(F)[ST]))) +
  theme_bw(base_size = 12)

figure4 <- p_fst_heat | p_ibd +
  plot_annotation(tag_levels = "A")

figure4

ggsave(file.path(figures_dir, "Figure4_FST_IBD.pdf"), figure4,
       width = 14, height = 6, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure4_FST_IBD.png"), figure4,
       width = 14, height = 6, dpi = 300)
