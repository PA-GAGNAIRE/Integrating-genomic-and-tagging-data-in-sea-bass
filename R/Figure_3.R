## FIGURE 3 ##

library('dplyr')
library('marmap')
library('ggpubr')
library('plotrix')

setwd("../data")

# Load bathymetry on the northeast Atlantic
ATLF <- getNOAA.bathy(lon1=-15,lon2=12,lat1=35,lat2=55,resolution=2)
# Compute transition object with path impossible in waters deeper than -200m depth
trans <- trans.mat(ATLF,min.depth=50,max.depth=-200)

# Load sample summer and winter coordinates (07/09/22)
COORDS <- read.table("SAMPLE_SUMMER_WINTER_COORDS.txt", header=T)
# Sort individuals by sample name 
COORDS <- arrange(COORDS, sample_filename)
# Load sample admixture proportions
ADMIX <- read.table("SAMPLE_ADMIXTURE_PROPORTIONS.txt", header=T)

# FIG 3A # Mediterranean ancestry gradient in First Summer

# Extract the SUMMER1 Coordinates, excluding NAs
SUMMER1_COORDS <- COORDS %>% dplyr::select(sample_name,Summer1_X,Summer1_Y) %>% drop_na(Summer1_X,Summer1_Y)

# Get the depth for each point, and filter by depth to remove points outside of the [50;-200] depth range
SUMMER1_DEPTH <- get.depth(ATLF, x=SUMMER1_COORDS$Summer1_X, y=SUMMER1_COORDS$Summer1_Y, locator=FALSE)
names(SUMMER1_DEPTH)[1] <- "Summer1_X"
names(SUMMER1_DEPTH)[2] <- "Summer1_Y"
SUMMER1_COORDS_DEPTH <- left_join(SUMMER1_COORDS, SUMMER1_DEPTH, by = c("Summer1_X" = "Summer1_X", "Summer1_Y" = "Summer1_Y"))
SUMMER1_COORDS <- SUMMER1_COORDS_DEPTH %>% filter(depth>-200) %>% filter(depth<50)

# Compute least cost distances for the transition matrix
ind_dist <- lc.dist(trans, rbind(c(-7.7991639, 37.00795),SUMMER1_COORDS[,2:3]), res ="dist")
SUMMER1_COORDS$dist <- ind_dist[1:nrow(SUMMER1_COORDS)]
names(SUMMER1_COORDS)[1] <- "IND"

# Join Summer1 individuals distances with sample admixture proportions
SUMMER1_ADMIX_DIST <- left_join(SUMMER1_COORDS, ADMIX, by = c("IND" = "IND"))

# Bin the distance values into quantiles (equal number of points per bin), using 20 bins
n_bins <- 20
bins <- quantile(SUMMER1_ADMIX_DIST$dist, probs = seq(0, 1, length.out = n_bins + 1))
x_bins <- cut(SUMMER1_ADMIX_DIST$dist, breaks = bins, include.lowest = TRUE)
# Calculate midpoints for each bin
midpoints <- sapply(1:(length(bins) - 1), function(i) mean(bins[i:(i+1)]))
# Map midpoints to each bin
mid_x <- midpoints[as.numeric(x_bins)]
# Create a list of y values corresponding to each bin
y_list <- split(SUMMER1_ADMIX_DIST$Q1, x_bins)

# Plot the boxplots using midpoints as x-axis values
myPal <- colorRampPalette(c("red","gold","blue"))
boxplot(y_list, at = midpoints, names = round(midpoints, 2), outline = F, xaxt = 'n', xlim=c(1320,2750), ylim=c(0,0.3), boxwex = 15, col=myPal(n_bins), xlab = "Interval midpoint distance to FA (km)", ylab = "Mediterranean ancestry")
rect(min(SUMMER1_ADMIX_DIST$dist),-0.05,1997,0.35,col = rgb(1, 0.5, 0, alpha = 0.1), border = NA)
rect(1997,-0.05,max(SUMMER1_ADMIX_DIST$dist)+100,0.35,col = rgb(0, 1, 0, alpha = 0.1), border = NA)
axis(1, at = pretty(SUMMER1_ADMIX_DIST$dist), labels = pretty(SUMMER1_ADMIX_DIST$dist))
abline(v=1729,lty=2,col="orange")
abline(v=2304,lty=2,col="green")
abline(v=2626,lty=2,col="green")
boxplot(y_list, at = midpoints, names = round(midpoints, 2), outline = F, xaxt = 'n', yaxt = 'n', boxwex = 15, col=myPal(n_bins), add = TRUE)
text(2000,0.29,expression(bold("First Summer")),cex=1.2)
text(1725,0.27,"BOB ICES stock area", col="orange")
text(1500,0.25,"27.8.b", col="orange")
text(1875,0.25,"27.8.a", col="orange")
text(2410,0.27,"NS ICES stock area", col="green")
text(2150,0.25,"27.7.e", col="green")
text(2150,0.235,"27.7.h", col="green")
text(2475,0.25,"27.7.d", col="green")
text(2725,0.25,"27.4.c", col="green")


# FIG 3C # Mediterranean ancestry gradient in First Winter

# Extract the WINTER1 Coordinates, excluding NAs
WINTER1_COORDS <- COORDS %>% dplyr::select(sample_name,Winter1_X,Winter1_Y) %>% drop_na(Winter1_X,Winter1_Y)

# Get the depth for each point, and filter by depth to remove points outside of the [50;-200] depth range
WINTER1_DEPTH <- get.depth(ATLF, x=WINTER1_COORDS$Winter1_X, y=WINTER1_COORDS$Winter1_Y, locator=FALSE)
names(WINTER1_DEPTH)[1] <- "Winter1_X"
names(WINTER1_DEPTH)[2] <- "Winter1_Y"
WINTER1_COORDS_DEPTH <- left_join(WINTER1_COORDS, WINTER1_DEPTH, by = c("Winter1_X" = "Winter1_X", "Winter1_Y" = "Winter1_Y"))
WINTER1_COORDS <- WINTER1_COORDS_DEPTH %>% filter(depth>-200) %>% filter(depth<50)

# Compute least cost distances for the transition matrix
ind_dist <- lc.dist(trans, rbind(c(-7.7991639, 37.00795),WINTER1_COORDS[,2:3]), res ="dist")
WINTER1_COORDS$dist <- ind_dist[1:nrow(WINTER1_COORDS)]
names(WINTER1_COORDS)[1] <- "IND"

# Join Winter1 individuals distances with sample admixture proportions
WINTER1_ADMIX_DIST <- left_join(WINTER1_COORDS, ADMIX, by = c("IND" = "IND"))

# Bin the distance values into quantiles (equal number of points per bin), using 20 bins
n_bins <- 20
bins <- quantile(WINTER1_ADMIX_DIST$dist, probs = seq(0, 1, length.out = n_bins + 1))
x_bins <- cut(WINTER1_ADMIX_DIST$dist, breaks = bins, include.lowest = TRUE)
# Calculate midpoints for each bin
midpoints <- sapply(1:(length(bins) - 1), function(i) mean(bins[i:(i+1)]))
# Map midpoints to each bin
mid_x <- midpoints[as.numeric(x_bins)]
# Create a list of y values corresponding to each bin
y_list <- split(WINTER1_ADMIX_DIST$Q1, x_bins)

# Plot the boxplots using midpoints as x-axis values
myPal <- colorRampPalette(c("red","gold","blue"))
boxplot(y_list, at = midpoints, names = round(midpoints, 2), outline = F, xaxt = 'n', xlim=c(1320,2750), ylim=c(0,0.3), boxwex = 15, col=myPal(n_bins), xlab = "Interval midpoint distance to FA (km)", ylab = "Mediterranean ancestry")
rect(min(WINTER1_ADMIX_DIST$dist),-0.05,1997,0.35,col = rgb(1, 0.5, 0, alpha = 0.1), border = NA)
rect(1997,-0.05,max(WINTER1_ADMIX_DIST$dist)+100,0.35,col = rgb(0, 1, 0, alpha = 0.1), border = NA)
axis(1, at = pretty(WINTER1_ADMIX_DIST$dist), labels = pretty(WINTER1_ADMIX_DIST$dist))
abline(v=1729,lty=2,col="orange")
abline(v=2304,lty=2,col="green")
abline(v=2626,lty=2,col="green")
boxplot(y_list, at = midpoints, names = round(midpoints, 2), outline = F, xaxt = 'n', yaxt = 'n', boxwex = 15, col=myPal(n_bins), add = TRUE)
text(2000,0.29,expression(bold("First Winter")),cex=1.2)
text(1725,0.27,"BOB ICES stock area", col="orange")
text(1500,0.25,"27.8.b", col="orange")
text(1875,0.25,"27.8.a", col="orange")
text(2410,0.27,"NS ICES stock area", col="green")
text(2150,0.25,"27.7.e", col="green")
text(2150,0.235,"27.7.h", col="green")
text(2475,0.25,"27.7.d", col="green")
text(2725,0.25,"27.4.c", col="green")
