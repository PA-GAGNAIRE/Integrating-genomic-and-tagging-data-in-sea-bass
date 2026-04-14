################################################################################
# Figure 5 – Electronic tagging: displacement trajectories and seasonal movements
# Gagnaire et al. (2026) Integrating genomic and tagging data reveals
# spatio-temporal population structure in Northeast Atlantic European sea bass
#
# Extracted from: Gagnaire_et_al_2026_figures.Rmd  ## FIGURE 5 ##
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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(lubridate)
library(geosphere)

# [CHANGE 3] Configurable paths – edit these two lines to match your setup
data_dir    <- "data"      # folder that contains the input CSV / TXT files
figures_dir <- "figures"   # folder where output figures will be saved

# [CHANGE 2] Load data
tagging <- read.csv(file.path(data_dir, "tagging_data.csv"))

# Create figures directory if it does not already exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- UNCHANGED CODE FROM RMD BELOW ----------

# Load coastline
world  <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world[world$continent == "Europe" | world$admin %in%
                  c("Morocco", "Algeria", "Tunisia", "Libya"), ]

# Tag and release information
release_sites <- tagging %>%
  distinct(tag_id, release_lon, release_lat, release_region)

# Recapture / detection events
detections <- tagging %>%
  filter(!is.na(recapture_lon)) %>%
  mutate(season = case_when(
    lubridate::month(as.Date(date)) %in% c(12, 1, 2) ~ "Winter",
    lubridate::month(as.Date(date)) %in% c(3, 4, 5)  ~ "Spring",
    lubridate::month(as.Date(date)) %in% c(6, 7, 8)  ~ "Summer",
    TRUE                                               ~ "Autumn"
  ))

# Displacement distance (Haversine)
detections <- detections %>%
  mutate(displacement_km = geosphere::distHaversine(
    cbind(release_lon, release_lat),
    cbind(recapture_lon, recapture_lat)
  ) / 1000)

season_colors <- c(Winter = "#4575B4", Spring = "#91CF60",
                   Summer = "#FC8D59", Autumn = "#D73027")

# --- Map of trajectories ---
p5_map <- ggplot() +
  geom_sf(data = europe, fill = "grey85", color = "grey60", linewidth = 0.3) +
  geom_segment(data = detections,
               aes(x = release_lon, y = release_lat,
                   xend = recapture_lon, yend = recapture_lat,
                   color = season),
               alpha = 0.4, linewidth = 0.5,
               arrow = arrow(length = unit(0.08, "inches"), type = "open")) +
  geom_point(data = release_sites,
             aes(x = release_lon, y = release_lat),
             shape = 21, fill = "yellow", color = "black",
             size = 3, stroke = 0.6) +
  scale_color_manual(values = season_colors, name = "Season") +
  coord_sf(xlim = c(-12, 10), ylim = c(38, 56), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "#D6EAF8"))

# --- Displacement distance by season ---
p5_disp <- ggplot(detections, aes(x = season, y = displacement_km, fill = season)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = season_colors) +
  labs(x = NULL, y = "Displacement distance (km)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

# --- Monthly displacement heatmap (individual × month) ---
monthly_disp <- detections %>%
  mutate(month = lubridate::month(as.Date(date), label = TRUE)) %>%
  group_by(tag_id, month) %>%
  summarise(mean_disp = mean(displacement_km, na.rm = TRUE), .groups = "drop")

p5_heat <- ggplot(monthly_disp, aes(x = month, y = tag_id, fill = mean_disp)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Mean disp. (km)", na.value = "grey90") +
  labs(x = "Month", y = "Individual tag") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

figure5 <- (p5_map | (p5_disp / p5_heat)) +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "A")

figure5

ggsave(file.path(figures_dir, "Figure5_tagging.pdf"), figure5,
       width = 14, height = 7, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Figure5_tagging.png"), figure5,
       width = 14, height = 7, dpi = 300)
