## FIGURE 2 ##

library('StAMPP')
library('adegenet')
library('parallel')
library('poppr')
library('hierfstat')
library('dendextend')
library('pheatmap')
library('viridis')
library('marmap')
library('tidyverse')
library('plotrix')

setwd("../data/")

# Load genotype data with the read.PLINK function to build an object of class genlight
gDAT <- read.PLINK(file='BARFRAY_765_ATL_SAMPLES__47211_SNPS_MAF0.01_HWE0.0001.raw', n.cores=1)

# Add chromosome and position information from plink .map file
gMAP <- read.table('BARFRAY_765_ATL_SAMPLES__47211_SNPS_MAF0.01_HWE0.0001.map', header=F)
chr(gDAT) <- gMAP$V1
position(gDAT) <- gMAP$V4

# Add tagging site information as population
gPOP <- read.table('BARFRAY_765_ATL_SAMPLES_POPMAP.txt', header=T)
pop(gDAT) <- gPOP[order(gPOP$IID),]$POP_CODE

# Print the sample size for each pop
summary(gDAT$pop)

# FIG 2A # Heatmap of pairwise Fst values with hierarchical clustering

# Compute pairwise FST (Weir & Cockerham 1984), with 10000 bootstraps across loci to generate 95% CIs and p-values
# pw_fst <- stamppFst(gDAT, nboots=10000, percent=95, nclusters=10)

# load Fst values and make a square matrix
fst <- pw_fst$Fsts
diag(fst) <- 0
symm_fst <- fst
symm_fst[upper.tri(symm_fst)] <- t(fst)[upper.tri(symm_fst)]
symm_fst <- pmax(symm_fst, 0)

# Make hierarchical clustering of pairwise Fst values
hclust_loc <- hclust(dist(symm_fst), method = "complete")

# Flip nodes to order samples locations (does not affect the tree)
loc_order <- c("FA","CB","OL","NO","LT","AD","SQ","SM","CH","PB","DK")
hclust_loc_col <- rotate(hclust_loc, loc_order)
hclust_loc_row <- rotate(hclust_loc, loc_order)

# Prepare annotation colors
tree_col <- data.frame(cutree(hclust_loc, k = 3))
colnames(tree_col) <- "cluster"
tree_col$cluster <- ifelse(tree_col$cluster==1, "cluster2",
                   ifelse(tree_col$cluster==2, "cluster3", "cluster1"))

ices_col <- as.factor(c("NS","NS","BOB","BOB","BOB","NS","BOB","NS","BOB","NS","Iberian"))
tree_col$ices <- ices_col

loc_col <- data.frame(location = loc_order)
row.names(loc_col) <- loc_col$location

myPal <- colorRampPalette(c("blue","gold","red"))
fill=c(myPal(10),"red4")

my_colour = list(
    location = c(FA=fill[11],CB=fill[10],OL=fill[9],NO=fill[8],LT=fill[7],AD=fill[6],SQ=fill[5],SM=fill[4],CH=fill[3],PB=fill[2],DK=fill[1]),
    ices = c(Iberian="purple", BOB="orange", NS="green"),
    cluster = c(cluster1="violet", cluster2="#ffb366", cluster3="lightgreen"))

# load P-values and make a matrix
Pv <- pw_fst$Pvalues
diag(Pv) <- 1
symm_Pv <- Pv
symm_Pv[upper.tri(symm_Pv)] <- t(Pv)[upper.tri(symm_Pv)]

# Function to assign significance levels
signif <- function(p)
  {ifelse(p == 1, ".",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", "NS"))))}

# Apply the function to the matrix of P-values
signif_matrix <- matrix(signif(symm_Pv), nrow = nrow(symm_Pv), ncol = ncol(symm_Pv))

# Make the plot
Q1 <- pheatmap(symm_fst, color = viridis(n = 256, alpha = 1, begin = 0.15, end = 1, direction = -1), breaks = seq(0, 0.001, length.out = 257), display_numbers = signif_matrix, fontsize_number = 10, cluster_rows = hclust_loc_row, cluster_cols = hclust_loc_col, annotation_colors = my_colour, annotation_row = tree_col, annotation_col = loc_col, treeheight_row = 0, treeheight_col = 150)


# FIG 2B # PCA analysis

# Perform a PCA analysis that includes all individuals (from the 10 French location and southern Portugal)
pca <- glPca(gDAT, nf=2, parallel=TRUE, n.cores=20)

# Compute variance explained by the 10 first PC axes
axis_var <- data.frame(axis = c(1:10), eig = pca$eig[c(1:10)]) %>% dplyr::mutate(p.eig = 100*(eig/sum(pca$eig)))
pc1 <- paste("PC1 (", paste(round(axis_var$p.eig[1],2)), "%)",sep="")
pc2 <- paste("PC2 (", paste(round(axis_var$p.eig[2],2)), "%)",sep="")

# Plot the PCA representing only the 10 populations from BARFRAY, with PC1 oriented by Mediterranean ancestry
coords <- pca$scores[1:705,]
coords[,1] <- -coords[,1]

myPal <- colorRampPalette(c("blue","gold","red"))
Q2 = with(pca, plot(coords, yaxt="n", xlim=c(-8,6), ylim=c(-12,10), ylab=pc2, xlab=pc1, col="white"))
legend("topleft",legend=c("CB","OL","NO","LT","AD","SQ","SM","CH","PB","DK") ,horiz=T,cex=1.05, fill=rev(myPal(10)))
axis(2, las=2)
abline(h=0,v=0,col="grey", lty=2)

s.class(coords, fac=factor(gPOP[order(gPOP$IID),][1:705,]$POP_CODE), 
        add.plot=T, 
        cpoint=2,
        clabel=0,
        pch=20,
        axesell=F,
        addaxes=F,
        cstar=0,
        col=transp(myPal(10)[c(6,10,3,1,7,9,8,2,4,5)],0.6),
        cellipse = 1.5)


# FIG 2C # Mediterranean ancestry gradient

# Load sample release coordinates
RELEAS <- read.table("SAMPLE_RELEASE_COORDS.txt", header=T)
POPX <- aggregate(x=RELEAS$X,by=list(RELEAS$population),FUN = mean)
colnames(POPX) <- c("POP","X")
POPY <- aggregate(x=RELEAS$Y,by=list(RELEAS$population),FUN = mean)
colnames(POPY) <- c("POP","Y")
POP_COORDS <- left_join(POPX,POPY,by='POP')
POP_COORDS <- POP_COORDS[c(2,7,6,5,1,10,9,3,8,4),]

# Load sample admixture proportions, calculate mean and SE per location
ADMIX <- read.table("SAMPLE_ADMIXTURE_PROPORTIONS.txt", header=T)
MEAN_MED <- aggregate(x=ADMIX$Q1,by=list(ADMIX$Pop),FUN = mean)
colnames(MEAN_MED) <- c("POP","MED_ANC")
MEAN_MED <- MEAN_MED[c(2,8,7,6,1,11,10,3,9,4),]
SE_MED <- aggregate(x=ADMIX$Q1,by=list(ADMIX$Pop),FUN = std.error)
colnames(SE_MED) <- c("POP","MED_ANC_SE")
SE_MED <- SE_MED[c(2,8,7,6,1,11,10,3,9,4),]

# Load bathymetry on the northeast Atlantic
ATLF <- getNOAA.bathy(lon1=-15,lon2=12,lat1=35,lat2=55,resolution=2)

# Compute transition object with path impossible in waters deeper than -200 meters depth
trans <- trans.mat(ATLF,min.depth=5,max.depth=-200)

# Compute least cost distances for the transition matrix
out_dist <- lc.dist(trans, rbind(c(-7.7991639, 37.00795),POP_COORDS[,2:3]), res ="dist")

# Add least cost distance to Faro location to the Mediterranean ancestry dataframe
MEAN_MED$dist <- out_dist[1:10]

# Make the plot
myPal <- colorRampPalette(c("red","gold","blue"))
plot(MED_ANC~dist, data=MEAN_MED, col=myPal(10), pch=20, cex=3, xlim=c(1500,2700),ylim=c(0,0.13), xlab="Distance to FA - southern Portugal (km)", ylab="Mean Mediterranean ancestry")

# Fit the tanh model
tanh_model <- nls(MED_ANC ~ a * tanh(b * rev(dist) + c) + d, data=MEAN_MED,
                  start = list(a = 1, b = 0.0005, c = -0.2, d = 0), algorithm = "default")

# Add the fitted curve for the observed data points
lines(MEAN_MED$dist, predict(tanh_model), col= "darkgrey", lty =2 , lwd = 2)
points(MED_ANC~dist, data=MEAN_MED, col=myPal(10), pch=20, cex=3)
segments(MEAN_MED$dist,MEAN_MED$MED_ANC-SE_MED$MED_ANC_SE,MEAN_MED$dist,MEAN_MED$MED_ANC+SE_MED$MED_ANC_SE, col=myPal(10), lwd = 2)
text(MEAN_MED$dist+30,MED_ANCESTRY+0.01,POP_COORDS$POP, col=myPal(10))
