################################################################################
# Figure 1 – Sampling map
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 1 ##
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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
metadata <- read.csv(file.path(data_dir, "metadata.csv"))

# Create figures directory if it does not already exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- UNCHANGED CODE FROM RMD BELOW ----------

# Define color palette for sampling regions
region_colors <- c(
  "Bay of Biscay"        = "#2166AC",
  "English Channel"      = "#4DAC26",
  "North Sea"            = "#1A9850",
  "Iberian Atlantic"     = "#D73027",
  "Western Mediterranean" = "#F46D43",
  "Eastern Mediterranean" = "#FDAE61"
)

# Load coastline and country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world[world$continent == "Europe" | world$admin %in%
                  c("Morocco", "Algeria", "Tunisia", "Libya"), ]

# Compute sample size per site for point scaling
site_summary <- metadata %>%
  group_by(site, longitude, latitude, region) %>%
  summarise(n = n(), .groups = "drop")

# Main map
fig1_map <- ggplot() +
  geom_sf(data = europe, fill = "grey85", color = "grey60", linewidth = 0.3) +
  geom_point(data = site_summary,
             aes(x = longitude, y = latitude,
                 fill = region, size = n),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.9) +
  geom_text_repel(data = site_summary,
                  aes(x = longitude, y = latitude, label = site),
                  size = 2.8, box.padding = 0.4, point.padding = 0.3,
                  segment.color = "grey50", segment.size = 0.3) +
  scale_fill_manual(values = region_colors, name = "Region") +
  scale_size_continuous(name = "Sample size", range = c(2, 8),
                        breaks = c(10, 30, 60, 100)) +
  coord_sf(xlim = c(-20, 40), ylim = c(28, 62), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "right",
    panel.background = element_rect(fill = "#D6EAF8"),
    panel.grid.major = element_line(color = "white", linewidth = 0.3)
  )

# Inset bar chart of sample sizes per region
fig1_inset <- ggplot(site_summary %>%
                       group_by(region) %>%
                       summarise(total = sum(n), .groups = "drop"),
                     aes(x = reorder(region, -total), y = total, fill = region)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = region_colors) +
  labs(x = NULL, y = "N individuals") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

# Combine map and inset
figure1 <- ggdraw(fig1_map) +
  draw_plot(fig1_inset, x = 0.60, y = 0.05, width = 0.38, height = 0.30)

figure1

ggsave(file.path(figures_dir, "Figure1_sampling_map.pdf"), figure1,
       width = 10, height = 7, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure1_sampling_map.png"), figure1,
       width = 10, height = 7, dpi = 300)
