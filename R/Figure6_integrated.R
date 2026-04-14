################################################################################
# Figure 6 – Integrated comparison of genomic ancestry and tagging connectivity
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 6 ##
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
library(ggrepel)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
metadata  <- read.csv(file.path(data_dir, "metadata.csv"))
admix_K2  <- read.table(file.path(data_dir, "admixture_K2.txt"))
tagging   <- read.csv(file.path(data_dir, "tagging_data.csv"))

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

# Attach metadata to admixture results (individuals must be in the same order)
admix_K2_df <- cbind(metadata, admix_K2)

# Recapture / detection events with displacement
detections <- tagging %>%
  filter(!is.na(recapture_lon)) %>%
  mutate(displacement_km = geosphere::distHaversine(
    cbind(release_lon, release_lat),
    cbind(recapture_lon, recapture_lat)
  ) / 1000)

# --- Atlantic ancestry proportion per site (from ADMIXTURE K = 2) ---
genomic_summary <- admix_K2_df %>%
  group_by(site, region, longitude) %>%
  summarise(
    mean_atlantic  = mean(V1),
    se_atlantic    = sd(V1) / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(longitude)

# --- Connectivity index from tagging (proportion of fish moving > 500 km) ---
tag_connectivity <- detections %>%
  group_by(release_region) %>%
  summarise(
    n_total    = n(),
    n_distant  = sum(displacement_km > 500, na.rm = TRUE),
    connect_idx = n_distant / n_total,
    se_connect  = sqrt(connect_idx * (1 - connect_idx) / n_total),
    .groups = "drop"
  )

# Merge genomic and tagging summaries on shared region
integrated <- genomic_summary %>%
  left_join(tag_connectivity, by = c("region" = "release_region"))

# --- Panel A: genomic ancestry by longitude ---
p6a <- ggplot(genomic_summary,
              aes(x = longitude, y = mean_atlantic, color = region)) +
  geom_errorbar(aes(ymin = mean_atlantic - se_atlantic,
                    ymax = mean_atlantic + se_atlantic),
                width = 0.5, alpha = 0.5) +
  geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(x = "Longitude (°E)", y = "Atlantic ancestry proportion") +
  theme_bw(base_size = 12)

# --- Panel B: tagging connectivity by longitude ---
p6b <- ggplot(integrated %>% filter(!is.na(connect_idx)),
              aes(x = longitude, y = connect_idx, color = region)) +
  geom_errorbar(aes(ymin = connect_idx - se_connect,
                    ymax = connect_idx + se_connect),
                width = 0.5, alpha = 0.5) +
  geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE,
              color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(x = "Longitude (°E)", y = "Tagging connectivity index") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# --- Panel C: scatter of genomic vs tagging across sites ---
p6c <- ggplot(integrated, aes(x = mean_atlantic, y = connect_idx, color = region)) +
  geom_errorbar(aes(ymin = connect_idx - se_connect,
                    ymax = connect_idx + se_connect),
                width = 0.01, alpha = 0.5) +
  geom_errorbarh(aes(xmin = mean_atlantic - se_atlantic,
                     xmax = mean_atlantic + se_atlantic),
                 height = 0.01, alpha = 0.5) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = site), size = 3, box.padding = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "black", fill = "grey80", linewidth = 0.8) +
  scale_color_manual(values = region_colors, name = "Region") +
  labs(x = "Atlantic ancestry proportion",
       y = "Tagging connectivity index") +
  theme_bw(base_size = 12)

figure6 <- (p6a / p6b | p6c) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

figure6

ggsave(file.path(figures_dir, "Figure6_integrated.pdf"), figure6,
       width = 14, height = 8, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure6_integrated.png"), figure6,
       width = 14, height = 8, dpi = 300)
