################################################################################
# Figure 2 – Principal Component Analysis (PCA) of genome-wide SNP data
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 2 ##
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
library(patchwork)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
geno     <- read.table(file.path(data_dir, "genotypes.txt"), header = TRUE, row.names = 1)
metadata <- read.csv(file.path(data_dir, "metadata.csv"))

# Create figures directory if it does not already exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- UNCHANGED CODE FROM RMD BELOW ----------

# Define color palette for sampling regions (shared across figures)
region_colors <- c(
  "Bay of Biscay"        = "#2166AC",
  "English Channel"      = "#4DAC26",
  "North Sea"            = "#1A9850",
  "Iberian Atlantic"     = "#D73027",
  "Western Mediterranean" = "#F46D43",
  "Eastern Mediterranean" = "#FDAE61"
)

# Perform PCA on the genotype matrix (coded as allele dosage 0/1/2)
# Replace missing values with column means before PCA
geno_imp <- apply(geno, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
pca_res <- prcomp(geno_imp, scale. = TRUE)

# Proportion of variance explained
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

# Combine PCA scores with metadata
pca_df <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  PC3    = pca_res$x[, 3],
  PC4    = pca_res$x[, 4],
  ID     = rownames(pca_res$x)
) %>%
  left_join(metadata, by = c("ID" = "individual_id"))

# PC1 vs PC2
p_12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region, shape = region)) +
  geom_point(alpha = 0.75, size = 2.5) +
  stat_ellipse(aes(group = region), type = "t", level = 0.9,
               linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 7), name = "Region") +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# PC3 vs PC4
p_34 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = region, shape = region)) +
  geom_point(alpha = 0.75, size = 2.5) +
  stat_ellipse(aes(group = region), type = "t", level = 0.9,
               linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = region_colors, name = "Region") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 7), name = "Region") +
  labs(x = paste0("PC3 (", var_exp[3], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)")) +
  theme_bw(base_size = 12)

# Scree plot (variance explained per PC)
p_scree <- data.frame(PC = 1:20, VarExp = var_exp[1:20]) %>%
  ggplot(aes(x = PC, y = VarExp)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_line(color = "darkred", linewidth = 0.7) +
  geom_point(color = "darkred", size = 2) +
  labs(x = "Principal Component", y = "Variance explained (%)") +
  theme_classic(base_size = 10)

figure2 <- (p_12 | p_34) / p_scree +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(tag_levels = "A")

figure2

ggsave(file.path(figures_dir, "Figure2_PCA.pdf"), figure2,
       width = 12, height = 8, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure2_PCA.png"), figure2,
       width = 12, height = 8, dpi = 300)
