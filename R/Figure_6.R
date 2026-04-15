## FIGURE 6 ##

library('ggplot2')
library('dartR')
library('sf')
library('tidyverse')
library('dplyr')
library('ggOceanMaps')
library('ggnewscale')
library('ggrepel')
library('ggpubr')
library('ggspatial')
library('ggforce')
library('scatterpie')

setwd("../data/")

## Dataset preparation ##

# Load genotype data from Taylor et al spawning samples in northeastern Atlantic
# 15 geographic samples and 41022 SNPs
Tspawn <- gl.load(file = "./Taylor/bass_reduced_ices_spawning.gl")

# Load genotype data from Taylor et al feeding samples in northeastern Atlantic
# 20 geographic samples and 41022 SNPs
Tfeed <- gl.load(file = "./Taylor/bass_reduced_ices_feeding.gl")

# Load Axiom_DlabCHIP Annotation File
Annot <- read.csv("Axiom_DlabCHIP_AnnotationFile.v1.csv", header = T)

# Extract Probe Set IDs from Taylor et al. and convert to dataframe
TIDs <- data.frame(do.call("rbind", strsplit(locNames(Tspawn), "_", fixed = TRUE)))
TIDs <- cbind(locNames(Tspawn), TIDs)
names(TIDs) <- c("Taylor_locNames", "Probe_Set_ID", "Allele_Taylor")

# Left join to Axiom_DlabCHIP Annotation File
Annot_TIDs <- left_join(Annot, TIDs)

# Extract Probe Set IDs from gDAT and convert to dataframe
IDs <- cbind(sub("(.*?_.*?)_.*", "\\1", locNames(gDAT)), data.frame(do.call("rbind", strsplit(locNames(gDAT), "_", fixed = TRUE))))
IDs <- cbind(locNames(gDAT), IDs)
names(IDs) <- c("locNames", "cust_id", "chr", "pos", "Allele")

# Left join to Axiom_DlabCHIP Annotation File
Annot_TIDs_IDs <- left_join(Annot_TIDs, IDs)

# Filter to only retain shared SNP loci between datasets
shared_snps <- Annot_TIDs_IDs %>% filter(!is.na(Taylor_locNames), !is.na(locNames))

# Sorts table first by chromosome, then by ascending SNP position within each chromosome
chr_levels <- c("LG1A", "LG1B", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12",
                "LG13", "LG14", "LG15","LG16", "LG17", "LG18-21", "LG19", "LG20", "LG22-25", "LG24", "LGx", "UN")
shared_snps_sorted <- shared_snps %>% mutate(cust_chr = factor(cust_chr, levels = chr_levels), cust_pos = as.numeric(cust_pos)) %>% arrange(cust_chr, cust_pos)

# Extract the whitelist of shared loci using both naming schemes for subsetting:
whitelist <- shared_snps_sorted$locNames
Twhitelist <- shared_snps_sorted$Taylor_locNames

# Subset each genlight object, preserving the order of loci:
Subset_Tspawn <- match(Twhitelist, locNames(Tspawn))
gl_Tspawn_sub <- Tspawn[, Subset_Tspawn]
locNames(gl_Tspawn_sub) <- whitelist

Subset_Tfeed <- match(Twhitelist, locNames(Tfeed))
gl_Tfeed_sub <- Tfeed[, Subset_Tfeed]
locNames(gl_Tfeed_sub) <- whitelist

Subset_gDAT <- match(whitelist, locNames(gDAT))
gl_gDAT_sub <- gDAT[, Subset_gDAT]

# Merge gl_Tspawn_sub with gl_Tfeed_sub
gl_Taylor_merged <- rbind(gl_Tspawn_sub, gl_Tfeed_sub)

# Sanity check assessing the correlation of allele frequencies with Taylor filtered dataset 
af_Taylor <- glMean(gl_Taylor_merged)
af_gDAT  <- glMean(gl_gDAT_sub)
af_df <- data.frame(locus = names(af_Taylor), af_Taylor = af_Taylor, af_gDAT = af_gDAT, whitelist = whitelist, Twhitelist = Twhitelist)
plot(af_df$af_Taylor, af_df$af_gDAT, xlab = "af_Taylor", ylab = "af_gDAT", pch  = 16, cex  = 0.6)
abline(0, 1, col = "red", lwd = 2)
mod <- lm(af_df$af_Taylor~af_df$af_gDAT)
abline(mod, col="green")

# Finally, merge the two filtered datasets
gl_merged <- rbind(gl_Taylor_merged, gl_gDAT_sub)
# Add chromosome and position information to Optional content and remove unnecessary loc.all
chr(gl_merged) <- shared_snps_sorted$cust_chr
position(gl_merged) <- shared_snps_sorted$cust_pos
gl_merged@loc.all <- NULL

# Print the sample size for each pop
summary(gl_merged$pop)

# FIG 6A # PCA analysis of Northeast Atlantic European seabass

# Perform a PCA analysis that includes all individuals (from Taylor's dataset plus the 10 French locations and southern Portugal)
pca <- glPca(gl_merged, nf=2, parallel=TRUE, n.cores=20)

# Compute variance explained by the 10 first PC axes
axis_var <- data.frame(axis = c(1:10), eig = pca$eig[c(1:10)]) %>% dplyr::mutate(p.eig = 100*(eig/sum(pca$eig)))
pc1 <- paste("PC1 (", paste(round(axis_var$p.eig[1],2)), "%)",sep="")
pc2 <- paste("PC2 (", paste(round(axis_var$p.eig[2],2)), "%)",sep="")

# Extract first two PC scores:
scores <- as.data.frame(pca$scores)

# Annotate population samples
scores$pop <- pop(gl_merged)
scores$pop2 <- c(factor(rep("Taylor", 844)), pop(gl_merged)[845:1609])
scores$season <- c(factor(rep("Spawning", 286)), factor(rep("Feeding", 558)), pop(gl_merged)[845:1609])

# Plot the PCA representing in color the 10 populations from BARFRAY overlapped with the admixed sample FA and Taylor's dataset, with PC1 oriented by Mediterranean ancestry

myPal <- colorRampPalette(c("blue","gold","red"))

S1 = with(pca, plot(scores, yaxt="n", ylab=pc2, xlab=pc1, xlim=c(-8,5), ylim=c(-8,12), col="white"))
legend("topleft",legend=c("CB","OL","NO","LT","AD","SQ","SM","CH","PB","DK","Taylor") ,horiz=T,cex=0.7, fill=c(rev(myPal(10)),"grey50"))
axis(2, las=2)
abline(h=0,v=0,col="grey50", lty=2)

s.class(scores[scores$pop2!="FA",], fac=factor(scores[scores$pop2!="FA",]$pop2), 
        add.plot=T, 
        cpoint=2,
        clabel=0,
        pch=20,
        axesell=F,
        addaxes=F,
        cstar=0,
        col=transp(c("grey50", myPal(10)[c(6,10,3,1,7,8,9, 2,4,5)]),0.6),
        cellipse = 1.5)


# FIG 6B # European Seabass from around the UK belong to the NS spawning stock

# Prepare data frame of individual PC1 coordinates:
sampID <- rownames(scores)
PC1 <- data.frame(cbind(sampID,scores$PC1))
names(PC1) <- c("sampID", "pc1")

# Load sample summer and winter coordinates (07/09/22) and extract winter1 positions
COORDS <- read.table("SAMPLE_SUMMER_WINTER_COORDS.txt", header=T)
COORDS <- arrange(COORDS, sample_filename)
POPS = c("DK","PB","CH","SM","SQ","AD","LT","OL","NO","CB")

IND_COORDS = NULL
for (p in 1:length(POPS))
  {
  sampID <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$sample_filename
  X <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_X
  Y <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_Y
  POP_P <- cbind(sampID,POPS[p],X,Y)
  IND_COORDS <- rbind(IND_COORDS,POP_P)
  }
IND_COORDS = data.frame(IND_COORDS)
colnames(IND_COORDS)  <- c("sampID", "POP", "X", "Y")

# Left join individual PC1 coordinates and winter1 positions:
PC1_XY <- left_join(PC1,IND_COORDS,by='sampID') %>% drop_na(X,Y)
i <- c(2, 4, 5)  
PC1_XY[, i] <- apply(PC1_XY[, i], 2, function(x) as.numeric(x))

# Define spawning stocks Reference Baselines based on winter latitude.
# Individuals south of latitude 48 are assigned to BOB,
# north of 48 to NS. These serve as Reference Baselines.
PC1_XY_Stock <- PC1_XY %>% mutate(Stock = if_else(Y < 48, "BOB", "NS"))

# Function to perform Kernel density estimation with bootstrap confidence intervals.
# This function estimates the empirical density of a PC1 coordinates
# and quantifies uncertainty using bootstrap resampling.
density_ci <- function(x, grid, nboot = 1000, conf = 0.95) {
  boot_dens <- replicate(
    nboot,
    {
      xs <- sample(x, replace = TRUE)
      density(xs, from = min(grid), to = max(grid), n = length(grid))$y
    }
  )
  alpha <- (1 - conf)/2
  list(
    y   = rowMeans(boot_dens),
    low = apply(boot_dens, 1, quantile, probs = alpha),
    high= apply(boot_dens, 1, quantile, probs = 1 - alpha)
  )
}

# Compute Reference Baselines PC1 densities for each spawning stock:
pc1_grid <- seq(-5, 5, length.out = 500)

dens_BOB <- density_ci(PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "BOB"], pc1_grid)
dens_NS <- density_ci(PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "NS"], pc1_grid)

# Assemble Reference Baselines' densities in a dataframe
dens_df <- bind_rows(
  data.frame(pc1 = pc1_grid, density = dens_BOB$y, low = dens_BOB$low, high = dens_BOB$high, group = "BOB"),
  data.frame(pc1 = pc1_grid, density = dens_NS$y, low = dens_NS$low, high = dens_NS$high, group = "NS"))

# Plot Reference Baselines' densities with empirical confidence intervals
# The figure 6B shows PC1 empirical distributions for BOB (orange) and NS (green) Reference Baselines (assigned with tagging data); alongside empirical 95% confidence intervals.
# Empirical distribution of PC1 coordinates for samples from Taylor et al is show with grey histogram.

Taylor_pc1 <- PC1[1:844,] %>% mutate(pc1 = as.numeric(pc1))
Taylor_pc1$group <- "Taylor"

S2 <- ggplot(dens_df, aes(pc1, density, color = group, fill = group)) +
  # Density lines
  geom_line(size = 1) +
  # Confidence intervals
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.25, color = NA) +
  # Histogram of individuals from Taylor
  geom_histogram(data = Taylor_pc1, aes(x = pc1, y = after_stat(density), fill = group),
  inherit.aes = FALSE, bins = 30, color = "grey50", alpha = 0.4) +
  # Stock colours  
  scale_color_manual(values = c("BOB" = "orange", "NS" = "green")) +
  scale_fill_manual(values = c("BOB" = "orange", "NS" = "green", "Taylor" = "grey50")) +
  theme_classic() +
  labs(x = "PC1 coordinate", y = "Density", fill = "Stock") +
  guides(color = "none")


## Empirical Bayesian mixed-stock framework ## 

# Store Reference Baselines' as density() objects
# These will be used for Bayesian inference of mixture proportions without individual assignment.
d_BOB <- density(PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "BOB"])
d_NS <- density(PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "NS"])

# Interpolate empirical densities
# Returns density values at any PC1 coordinates.
# rule = 2 avoids zero-density artefacts at distribution tails.
density_at <- function(x, dens) {
  approx(dens$x, dens$y, xout = x, rule = 2)$y
}

# Empirical Bayesian inference of mixture proportions
# Estimates the proportion 'pi' of BOB stock in a mixture sample (proportion of NS stock is 1-pi).
# Individuals are treated as independent draws from a two-stock mixture.
# Returns the full posterior distribution of pi.
estimate_mixture_proportion <- function(
  pc1_mix,  # PC1 coordinates of the mixture sample
  d_BOB,
  d_NS,
  prior = c(1, 1),  # Beta prior probability of each group (i.e. equal prior probability to belong to either stock)
  grid_length = 1000,
  ci = 0.95
) {
  # Grid of sample mixture proportions
  pi_grid <- seq(0, 1, length.out = grid_length)

  # Log-likelihood for each pi
  loglik <- sapply(pi_grid, function(pi) {
    probs <- pi * density_at(pc1_mix, d_BOB) +
             (1 - pi) * density_at(pc1_mix, d_NS)
    sum(log(probs))
  })

  # Log-prior Beta to use a uniform prior on the log scale
  # the prior contributes no information 
  logprior <- dbeta(pi_grid, prior[1], prior[2], log = TRUE) 

  # Log-posterior
  # the posterior is driven entirely by the likelihood (Log-prior non-informative)
  logpost <- loglik + logprior # Bayes’ theorem taking logs
  post <- exp(logpost - max(logpost)) # Convert back from log scale and subtract the maximum to ensure the largest exponentiated value is 1
  post <- post / sum(post) # normalize to make a proper probability distribution

  # Posterior summaries
  alpha <- (1 - ci) / 2
  list(
    pi_map  = pi_grid[which.max(post)], # Maximum A Posteriori (MAP) estimate of mixture proportion (point estimate)
    ci = c(
      pi_grid[which(cumsum(post) >= alpha)[1]],
      pi_grid[which(cumsum(post) >= 1 - alpha)[1]]
    ),
    posterior = data.frame(pi = pi_grid, density = post)
  )
}

# Bootstrap Reference Baselines' densities to propagate baseline uncertainty
# Reference individuals are resampled with replacement.
# For each bootstrap replicate:
#   - Reference Baselines' densities are recomputed
#   - mixture proportion is re-estimated
# Resulting CIs incorporate reference-density uncertainty.
estimate_mixture_bootstrap <- function(
  pc1_mix,
  pc1_BOB,
  pc1_NS,
  B = 1000,
  prior = c(1, 1), # Beta prior probability of each group (i.e. equal prior probability to belong to either stock) 
  grid_length = 1000,
  ci = 0.95
) {
  pi_mean_boot <- numeric(B)

  for (b in seq_len(B)) {

    # Bootstrap Baseline Reference samples to recompute densities
    d_BOB_b <- density(sample(pc1_BOB, replace = TRUE))
    d_NS_b  <- density(sample(pc1_NS, replace = TRUE))

    # Estimate mixture proportion
    res <- estimate_mixture_proportion(
      pc1_mix = pc1_mix,
      d_BOB = d_BOB_b,
      d_NS  = d_NS_b,
      prior = prior,
      grid_length = grid_length,
      ci = ci
    )

    pi_mean_boot[b] <- res$pi_map # sampling distribution of pi obtained by repeatedly resampling the reference baseline
  }

  alpha <- (1 - ci) / 2
  list(
    pi_mean = mean(pi_mean_boot),
    ci = quantile(
      pi_mean_boot,
      probs = c(alpha, 1 - alpha) # confidence interval based on the sampling distribution of pi (uncertainty due to baseline estimation)
    ),
    bootstrap_estimates = pi_mean_boot
  )
}


## Empirical Bayesian mixed-stock Analysis ## 

# Apply to Taylor's dataset treated as a mixture sample
result_ref <- estimate_mixture_proportion(
  pc1_mix = Taylor_pc1$pc1,
  d_BOB = d_BOB,
  d_NS  = d_NS
)

# Bootstrap-propagated uncertainty
result_boot <- estimate_mixture_bootstrap(
  pc1_mix = Taylor_pc1$pc1,
  pc1_BOB = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "BOB"],
  pc1_NS  = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "NS"],
  B = 1000
)

# FIG 6C #  Mixed-stock analysis within ICES area in SUMMER1

# Extract summer1 positions
IND_COORDS = NULL
for (p in 1:length(POPS))
  {
  sampID <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Summer1_X,Summer1_Y))$sample_filename
  X <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Summer1_X,Summer1_Y))$Summer1_X
  Y <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Summer1_X,Summer1_Y))$Summer1_Y
  POP_P <- cbind(sampID,POPS[p],X,Y)
  IND_COORDS <- rbind(IND_COORDS,POP_P)
  }
IND_COORDS = data.frame(IND_COORDS)
colnames(IND_COORDS)  <- c("sampID", "POP", "X", "Y")

# Left join individual PC1 coordinates and summer1 positions:
SUMMER1_PC1_XY <- left_join(PC1,IND_COORDS,by='sampID') %>% drop_na(X,Y)
i <- c(2, 4, 5)  
SUMMER1_PC1_XY[, i] <- apply(SUMMER1_PC1_XY[, i], 2, function(x) as.numeric(x))

# load ICES polygons and rectangles
data("ices_areas", package = "ggOceanMaps")

# Convert individual summer1 positions dataframe to sf (WGS84)
SUMMER1_PC1_XY_sf <- st_as_sf(SUMMER1_PC1_XY, coords = c("X", "Y"), crs = 4326)
# Transform individual positions to match ICES CRS (3995)
SUMMER1_PC1_XY_sf <- st_transform(SUMMER1_PC1_XY_sf, st_crs(ices_areas))
# Spatial join
SUMMER1_PC1_XY_sf_ICES <- st_join(SUMMER1_PC1_XY_sf, ices_areas, join = st_intersects)
# Drop geometry column to convert to a regular data frame
SUMMER1_PC1_XY_sf_ICES <- SUMMER1_PC1_XY_sf_ICES %>% st_drop_geometry()
# List of summer1 ICES areas for the map
SUMMER1_ICES <- c("27.8.c", "27.8.b", "27.8.a", "27.7.e", "27.7.d", "27.4.c")
# Filter summer1 ICES areas and select columns
SUMMER1_PC1_ICES <- SUMMER1_PC1_XY_sf_ICES %>% filter(Area_Full %in% SUMMER1_ICES) %>% dplyr::select(sampID, pc1, POP, Area_Full)

# Extract Feeding samples and PC1 coordinates from Taylor's dataset
FEEDING_PC1_POP <- scores %>% filter(season=="Feeding") %>% dplyr::select(PC1, pop)
colnames(FEEDING_PC1_POP) <- c("pc1", "ICESNAME")
# Read ICES rectangle shapefile from local disk (requires .shp, .dbf and ..shx files in the same directory)
rects_sf <- st_read(dsn = "./Taylor/ICES_Statistical_Rectangles_Eco.shp")
# Assign CRS for the rectangles
st_crs(rects_sf) <- 4326
# Make CRS match between rectangles and ICES areas so that they can be spatially joined
rects_sf <- st_transform(rects_sf, st_crs(ices_areas))
# Assign rectangles to areas
rects_with_area <- st_join(rects_sf, ices_areas, join = st_within, largest = TRUE)
# Select columns of rectangle and ICES areas' names
rects_df <- rects_with_area %>% st_drop_geometry() %>% dplyr::select(ICESNAME, Area_Full)
# Assign ICES areas to Taylor's Feeding dataset
FEEDING_PC1_ICES <- left_join(FEEDING_PC1_POP, rects_df, by = "ICESNAME")
# List of FEEDING ICES areas for the map
FEEDING_ICES <- c("27.7.g", "27.7.f", "27.7.d", "27.7.b", "27.7.a", "27.6.a", "27.4.c", "27.4.b", "27.3.a.20")
# Filter FEEDING ICES areas and select columns
FEEDING_PC1_ICES <- FEEDING_PC1_ICES %>% filter(Area_Full %in% FEEDING_ICES) %>% dplyr::select(pc1, ICESNAME, Area_Full)

# Vector of all ICES areas covered by FEEDING samples (merged dataset)
MERGED_FEEDING_ICES <- c("27.8.c", "27.8.b", "27.8.a", "27.7.e", "27.7.g", "27.7.f", "27.7.d", "27.7.b", "27.7.a", "27.6.a", "27.4.c", "27.4.b", "27.3.a.20")

FEEDING_ICES_MIX <- data.frame(
  Area_Full = character(),
  n = numeric(),
  BOB_Pi = numeric(),
  BOB_Pi_low = numeric(),
  BOB_Pi_hi = numeric(),
  BOB_Pi_boot = numeric(),
  BOB_Pi_boot_low = numeric(),
  BOB_Pi_boot_hi = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(MERGED_FEEDING_ICES)) {

  area_i <- MERGED_FEEDING_ICES[i]

  PC1_ICES <- bind_rows(
    SUMMER1_PC1_ICES %>% filter(Area_Full == area_i) %>% dplyr::select(pc1),
    FEEDING_PC1_ICES %>% filter(Area_Full == area_i) %>% dplyr::select(pc1)
  )
  
  Pi_ICES <- estimate_mixture_proportion(pc1_mix = PC1_ICES$pc1, d_BOB = d_BOB, d_NS  = d_NS)
  
  Pi_ICES_boot <- estimate_mixture_bootstrap(pc1_mix = PC1_ICES$pc1, 
                                             pc1_BOB = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "BOB"], 
                                             pc1_NS  = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "NS"], 
                                             B = 1000)
  
  FEEDING_ICES_MIX <- bind_rows(
    FEEDING_ICES_MIX,
    data.frame(
      Area_Full = area_i,
      n = nrow(PC1_ICES),
      BOB_Pi = Pi_ICES$pi_map,
      BOB_Pi_low = Pi_ICES$ci[1],
      BOB_Pi_hi = Pi_ICES$ci[2],
      BOB_Pi_boot = Pi_ICES_boot$pi_mean,
      BOB_Pi_boot_low = Pi_ICES_boot$ci[1],
      BOB_Pi_boot_hi = Pi_ICES_boot$ci[2]    
    )
  )
}

# Compute NS proportions
FEEDING_ICES_MIX <- FEEDING_ICES_MIX %>% mutate(
                                  # Compute NS proportions
                                  NS_Pi = 1-BOB_Pi, 
                                  NS_Pi_low = 1-BOB_Pi_low, 
                                  NS_Pi_hi = 1-BOB_Pi_hi,
                                  NS_Pi_boot = 1-BOB_Pi_boot, 
                                  NS_Pi_boot_low = 1-BOB_Pi_boot_low, 
                                  NS_Pi_boot_hi = 1-BOB_Pi_boot_hi, 
                                  # Pie radius scaled by sample size
                                  radius = log(n + 1) * 15000,
                                  # convert CI to radians
                                  arc_start = 2*pi*BOB_Pi_low,
                                  arc_end   = 2*pi*BOB_Pi_hi,
                                  # CI color depending on predominant origin
                                  CI_col = ifelse(BOB_Pi >= 0.5, "BOB", "NS"))

# Compute ICES centroids for positioning pies
ices_clean <- sf::st_transform(sf::st_make_valid(ices_areas), 4326)

centroids <- ices_clean %>%
  filter(Area_Full %in% FEEDING_ICES_MIX$Area_Full) %>%
  st_centroid() %>%
  cbind(sf::st_coordinates(.)) %>%
  st_drop_geometry() %>%
  select(Area_Full, X = X, Y = Y)

# Join with inferred mixture proportions dataframe:
FEEDING_ICES_MIX <- FEEDING_ICES_MIX %>% left_join(centroids, by = "Area_Full")

# Make coordinates match the map projection:
FEEDING_ICES_MIX_sf <- st_as_sf(FEEDING_ICES_MIX, coords = c("X","Y"), crs = 4326)
FEEDING_ICES_MIX_sf <- st_transform(FEEDING_ICES_MIX_sf, 3995)

# Extract the projected coordinates for plotting:
coords <- st_coordinates(FEEDING_ICES_MIX_sf)
FEEDING_ICES_MIX$x <- coords[,1]
FEEDING_ICES_MIX$y <- coords[,2]

# Make the map and overlay pie charts:
S3 = basemap(limits = c(-12,10,42,58), bathymetry = TRUE,  land.col = "#eeeac4", grid.col = NA, grid.size = 0.5, base_size = 18, bathy.style = "poly_greys") +
  labs(x = "Longitude", y = "Latitude") +
  ggnewscale::new_scale_fill() +
  annotation_spatial(ices_areas[36:37,], mapping = aes(), fill="purple", col="purple", alpha = 0.1) + 
  annotation_spatial(ices_areas[39:40,], mapping = aes(), fill="orange", col="orange", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(7:8,43:48),], mapping = aes(), fill="green", col="green", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(27,41:42),], mapping = aes(), fill="blue", col="blue", alpha = 0.1) + 
  theme(legend.position = "none")

S3 = reorder_layers(S3) + annotation_scale(location = "br") +
ggtitle("Mixture Proportions During Feeding Season")

S3 <- S3 + geom_scatterpie(data = FEEDING_ICES_MIX,
    aes(x = x, y = y, r = radius),
    cols = c("BOB_Pi","NS_Pi"),
    color = "black") +
  scale_fill_manual(values = c(
    BOB_Pi = "orange",
    NS_Pi  = "green")) +
  geom_arc(data = FEEDING_ICES_MIX,
    aes(x0 = x, y0 = y, r = radius * 1.15, start = arc_start, end = arc_end, colour = CI_col),
    linewidth = 1.4, alpha = 0.9) + 
  scale_colour_manual(
    values = c(BOB="green", NS="orange"),
    guide = "none")


# FIG 6D #  Mixed-stock analysis within ICES area in WINTER1

# Extract winter1 positions
IND_COORDS = NULL
for (p in 1:length(POPS))
  {
  sampID <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$sample_filename
  X <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_X
  Y <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_Y
  POP_P <- cbind(sampID,POPS[p],X,Y)
  IND_COORDS <- rbind(IND_COORDS,POP_P)
  }
IND_COORDS = data.frame(IND_COORDS)
colnames(IND_COORDS)  <- c("sampID", "POP", "X", "Y")

# Left join individual PC1 coordinates and winter1 positions:
WINTER1_PC1_XY <- left_join(PC1,IND_COORDS,by='sampID') %>% drop_na(X,Y)
i <- c(2, 4, 5)  
WINTER1_PC1_XY[, i] <- apply(WINTER1_PC1_XY[, i], 2, function(x) as.numeric(x))

# load ICES polygons and rectangles
data("ices_areas", package = "ggOceanMaps")

# Convert individual winter1 positions dataframe to sf (WGS84)
WINTER1_PC1_XY_sf <- st_as_sf(WINTER1_PC1_XY, coords = c("X", "Y"), crs = 4326)
# Transform individual positions to match ICES CRS (3995)
WINTER1_PC1_XY_sf <- st_transform(WINTER1_PC1_XY_sf, st_crs(ices_areas))
# Spatial join
WINTER1_PC1_XY_sf_ICES <- st_join(WINTER1_PC1_XY_sf, ices_areas, join = st_intersects)
# Drop geometry column to convert to a regular data frame
WINTER1_PC1_XY_sf_ICES <- WINTER1_PC1_XY_sf_ICES %>% st_drop_geometry()
# List of winter1 ICES areas for the map
WINTER1_ICES <- c("27.8.c", "27.8.b", "27.8.a", "27.7.h", "27.7.e", "27.7.d", "27.4.c")
# Filter winter1 ICES areas and select columns
WINTER1_PC1_ICES <- WINTER1_PC1_XY_sf_ICES %>% filter(Area_Full %in% WINTER1_ICES) %>% dplyr::select(sampID, pc1, POP, Area_Full)

# Extract Spawning samples and PC1 coordinates from Taylor's dataset
SPAWNING_PC1_POP <- scores %>% filter(season=="Spawning") %>% dplyr::select(PC1, pop)
colnames(SPAWNING_PC1_POP) <- c("pc1", "ICESNAME")
# Read ICES rectangle shapefile from local disk (requires .shp, .dbf and ..shx files in the same directory)
rects_sf <- st_read(dsn = "./Taylor/ICES_Statistical_Rectangles_Eco.shp")
# Assign CRS for the rectangles
st_crs(rects_sf) <- 4326
# Make CRS match between rectangles and ICES areas so that they can be spatially joined
rects_sf <- st_transform(rects_sf, st_crs(ices_areas))
# Assign rectangles to areas
rects_with_area <- st_join(rects_sf, ices_areas, join = st_within, largest = TRUE)
# Select columns of rectangle and ICES areas' names
rects_df <- rects_with_area %>% st_drop_geometry() %>% dplyr::select(ICESNAME, Area_Full)
# Assign ICES areas to Taylor's Feeding dataset
SPAWNING_PC1_ICES <- left_join(SPAWNING_PC1_POP, rects_df, by = "ICESNAME")
# List of FEEDING ICES areas for the map
SPAWNING_ICES <- c("27.7.f", "27.7.e", "27.7.d", "27.7.b", "27.7.a", "27.4.c")
# Filter FEEDING ICES areas and select columns
SPAWNING_PC1_ICES <- SPAWNING_PC1_ICES %>% filter(Area_Full %in% SPAWNING_ICES) %>% dplyr::select(pc1, ICESNAME, Area_Full)


# Vector of all ICES areas covered by SPAWNING samples (merged dataset)
MERGED_SPAWNING_ICES <- c("27.8.c", "27.8.b", "27.8.a", "27.7.e", "27.7.f", "27.7.d", "27.7.b", "27.7.a", "27.4.c")

SPAWNING_ICES_MIX <- data.frame(
  Area_Full = character(),
  n = numeric(),
  BOB_Pi = numeric(),
  BOB_Pi_low = numeric(),
  BOB_Pi_hi = numeric(),
  BOB_Pi_boot = numeric(),
  BOB_Pi_boot_low = numeric(),
  BOB_Pi_boot_hi = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(MERGED_SPAWNING_ICES)) {

  area_i <- MERGED_SPAWNING_ICES[i]

  PC1_ICES <- bind_rows(
    WINTER1_PC1_ICES %>% filter(Area_Full == area_i) %>% dplyr::select(pc1),
    SPAWNING_PC1_ICES %>% filter(Area_Full == area_i) %>% dplyr::select(pc1)
  )
  
  Pi_ICES <- estimate_mixture_proportion(pc1_mix = PC1_ICES$pc1, d_BOB = d_BOB, d_NS  = d_NS)
  
  Pi_ICES_boot <- estimate_mixture_bootstrap(pc1_mix = PC1_ICES$pc1, 
                                             pc1_BOB = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "BOB"], 
                                             pc1_NS  = PC1_XY_Stock$pc1[PC1_XY_Stock$Stock == "NS"], 
                                             B = 1000)
  
  SPAWNING_ICES_MIX <- bind_rows(
    SPAWNING_ICES_MIX,
    data.frame(
      Area_Full = area_i,
      n = nrow(PC1_ICES),
      BOB_Pi = Pi_ICES$pi_map,
      BOB_Pi_low = Pi_ICES$ci[1],
      BOB_Pi_hi = Pi_ICES$ci[2],
      BOB_Pi_boot = Pi_ICES_boot$pi_mean,
      BOB_Pi_boot_low = Pi_ICES_boot$ci[1],
      BOB_Pi_boot_hi = Pi_ICES_boot$ci[2]    
    )
  )
}

# Compute NS proportions
SPAWNING_ICES_MIX <- SPAWNING_ICES_MIX %>% mutate(
                                  # Compute NS proportions
                                  NS_Pi = 1-BOB_Pi, 
                                  NS_Pi_low = 1-BOB_Pi_low, 
                                  NS_Pi_hi = 1-BOB_Pi_hi,
                                  NS_Pi_boot = 1-BOB_Pi_boot, 
                                  NS_Pi_boot_low = 1-BOB_Pi_boot_low, 
                                  NS_Pi_boot_hi = 1-BOB_Pi_boot_hi, 
                                  # Pie radius scaled by sample size
                                  radius = log(n + 1) * 15000,
                                  # convert CI to radians
                                  arc_start = 2*pi*BOB_Pi_low,
                                  arc_end   = 2*pi*BOB_Pi_hi,
                                  # CI color depending on predominant origin
                                  CI_col = ifelse(BOB_Pi >= 0.5, "BOB", "NS"))

# Compute ICES centroids for positioning pies
ices_clean <- sf::st_transform(sf::st_make_valid(ices_areas), 4326)

centroids <- ices_clean %>%
  filter(Area_Full %in% SPAWNING_ICES_MIX$Area_Full) %>%
  st_centroid() %>%
  cbind(sf::st_coordinates(.)) %>%
  st_drop_geometry() %>%
  select(Area_Full, X = X, Y = Y)

# Join with inferred mixture proportions dataframe:
SPAWNING_ICES_MIX <- SPAWNING_ICES_MIX %>% left_join(centroids, by = "Area_Full")

# Make coordinates match the map projection:
SPAWNING_ICES_MIX_sf <- st_as_sf(SPAWNING_ICES_MIX, coords = c("X","Y"), crs = 4326)
SPAWNING_ICES_MIX_sf <- st_transform(SPAWNING_ICES_MIX_sf, 3995)

# Extract the projected coordinates for plotting:
coords <- st_coordinates(SPAWNING_ICES_MIX_sf)
SPAWNING_ICES_MIX$x <- coords[,1]
SPAWNING_ICES_MIX$y <- coords[,2]

# Make the map and overlay pie charts:
S4 = basemap(limits = c(-12,10,42,58), bathymetry = TRUE,  land.col = "#eeeac4", grid.col = NA, grid.size = 0.5, base_size = 18, bathy.style = "poly_greys") +
  labs(x = "Longitude", y = "Latitude") +
  ggnewscale::new_scale_fill() +
  annotation_spatial(ices_areas[36:37,], mapping = aes(), fill="purple", col="purple", alpha = 0.1) + 
  annotation_spatial(ices_areas[39:40,], mapping = aes(), fill="orange", col="orange", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(7:8,43:48),], mapping = aes(), fill="green", col="green", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(27,41:42),], mapping = aes(), fill="blue", col="blue", alpha = 0.1) + 
  theme(legend.position = "none")

S4 = reorder_layers(S4) + annotation_scale(location = "br") +
ggtitle("Mixture Proportions During Spawning Season")

S4 <- S4 + geom_scatterpie(data = SPAWNING_ICES_MIX,
    aes(x = x, y = y, r = radius),
    cols = c("BOB_Pi","NS_Pi"),
    color = "black") +
  scale_fill_manual(values = c(
    BOB_Pi = "orange",
    NS_Pi  = "green")) +
  geom_arc(data = SPAWNING_ICES_MIX,
    aes(x0 = x, y0 = y, r = radius * 1.15, start = arc_start, end = arc_end, colour = CI_col),
    linewidth = 1.4, alpha = 0.9) + 
  scale_colour_manual(
    values = c(BOB="green", NS="orange"),
    guide = "none")
