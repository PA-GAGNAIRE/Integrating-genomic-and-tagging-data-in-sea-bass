## FIGURE 4 ##

library('adegenet')
library('parallel')
library('tidyverse')
library('lme4')

setwd("../data")

# FIG 4 # Test for regional philopatry using a Generalized Linear Model (GLM) or a mixed-effects model (GLMM) to test for a significant correlation between individual coordinate on PC1 and Winter1 latitude across multiple tagging locations, while accounting for variation in the correlation across these locations. 

# Load genotype data with the read.PLINK function to build an object of class genlight
gDAT <- read.PLINK(file='BARFRAY_765_ATL_SAMPLES__47211_SNPS_MAF0.01_HWE0.0001.raw', n.cores=1)

# Add chromosome and position information from plink .map file
gMAP <- read.table('BARFRAY_765_ATL_SAMPLES__47211_SNPS_MAF0.01_HWE0.0001.map', header=F)
chr(gDAT) <- gMAP$V1
position(gDAT) <- gMAP$V4

# Add tagging site information as population
gPOP <- read.table('BARFRAY_765_ATL_SAMPLES_POPMAP.txt', header=T)
pop(gDAT) <- gPOP[order(gPOP$IID),]$POP_CODE

# Perform a PCA analysis that includes all individuals (from the 10 French location and southern Portugal)
pca <- glPca(gDAT, nf=2, parallel=TRUE, n.cores=20)

# Extract individual coordinates from FIG 2B 
coords <- pca$scores[1:705,]
pc1 <- -coords[,1]

# Prepare data frame of individual PC1 coordinates:
sampID <- names(pc1)
PC1 <- data.frame(cbind(sampID,pc1))

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

# look at the data split by tagging location
ggplot(aes(pc1, Y), data = PC1_XY) + 
  geom_point() + 
  facet_wrap(~ POP)

# Retain only NS tagging locations with more than one observation
PC1_XY_NS <- PC1_XY %>% filter(POP=="DK"|POP=="PB"|POP=="CH"|POP=="SQ")

# fit a random-slope and random-intercept mixed model to test whether genetic makeup (pc1) affects latitudinal spawning coordinate (Y), while assuming that tagging locations have different latitudinal spawning baselines and that the relationship with genetics may vary among tagging locations. 
mixed.ranslope <- lmer(Y ~ pc1 + (1 + pc1|POP), data = PC1_XY_NS)

# summarize the results
summary(mixed.ranslope)
stargazer(mixed.ranslope, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# Prepare the plot:
myPal <- colorRampPalette(c("blue","gold","red"))
fill_colors <- transp(myPal(10)[c(5,3,2,1)], alpha=0.5)

g <- ggplot(PC1_XY_NS, aes(x = pc1, y = Y, colour = POP)) +
      facet_wrap(~factor(POP, c("SQ","CH","PB","DK")),nrow=1) +
      geom_line(data = cbind(PC1_XY_NS, pred = predict(mixed.ranslope)), aes(y = pred), linewidth = 1) +
      scale_color_manual(values=c("SQ"=myPal(10)[5],"CH"=myPal(10)[3],"PB"=myPal(10)[2],"DK"=myPal(10)[1])) +
      geom_hline(yintercept=48, linetype=2) +
      geom_point(alpha=0.5, size=2.5, colour=myPal(10)[c(3,1,2,5)][as.numeric(factor(PC1_XY_NS$POP))]) +
      theme_classic() +
      labs(x = "PC1", y = "Spawning latitude") +
      theme(legend.position = "none", panel.spacing = unit(2, "lines"))

# Load sample migration phenotypes to add information on individual spawning site fidelity
PHENO = read.table("SAMPLE_MIGRATION_PHENOTYPE.txt", header=T)

FID_NS <- PHENO %>% na.exclude(spawning_ground_strategy) %>% filter(population=="DK"|population=="PB"|population=="CH"|population=="SQ") %>% select(sampID,spawning_ground_strategy)

# Left join spawning site fidelity with individual PC1 coordinates and winter1 positions:
FID_PC1_Y_NS <- left_join(FID_NS,PC1_XY_NS,by='sampID')
i <- c(3, 5, 6)  
FID_PC1_Y_NS[, i] <- apply(FID_PC1_Y_NS[, i], 2, function(x) as.numeric(x))

# Extract dataframes of individuals showing spawning site fidelity to each site for each tagging location
SQ_FID_BOB <- FID_PC1_Y_NS %>% filter(POP=="SQ") %>% filter(spawning_ground_strategy=="BOB") %>% select(pc1,Y,POP)
SQ_FID_NS <- FID_PC1_Y_NS %>% filter(POP=="SQ") %>% filter(spawning_ground_strategy=="ICCNS") %>% select(pc1,Y,POP)
CH_FID_BOB <- FID_PC1_Y_NS %>% filter(POP=="CH") %>% filter(spawning_ground_strategy=="BOB") %>% select(pc1,Y,POP)
CH_FID_NS <- FID_PC1_Y_NS %>% filter(POP=="CH") %>% filter(spawning_ground_strategy=="ICCNS") %>% select(pc1,Y,POP)
PB_FID_BOB <- FID_PC1_Y_NS %>% filter(POP=="PB") %>% filter(spawning_ground_strategy=="BOB") %>% select(pc1,Y,POP)
PB_FID_NS <- FID_PC1_Y_NS %>% filter(POP=="PB") %>% filter(spawning_ground_strategy=="ICCNS") %>% select(pc1,Y,POP)
DK_FID_BOB <- FID_PC1_Y_NS %>% filter(POP=="DK") %>% filter(spawning_ground_strategy=="BOB") %>% select(pc1,Y,POP)
DK_FID_NS <- FID_PC1_Y_NS %>% filter(POP=="DK") %>% filter(spawning_ground_strategy=="ICCNS") %>% select(pc1,Y,POP)

# Add individuals showing spawning site fidelity to the plot 
g <- g + geom_point(data = SQ_FID_BOB, size=2.5, shape=21, colour="orange") +
  geom_point(data = SQ_FID_NS, size=2.5, shape=21, colour="green") +
  geom_point(data = CH_FID_BOB, size=2.5, shape=21, colour="orange") +
  geom_point(data = CH_FID_NS, size=2.5, shape=21, colour="green") +
  geom_point(data = PB_FID_BOB, size=2.5, shape=21, colour="orange") +
  geom_point(data = PB_FID_NS, size=2.5, shape=21, colour="green") +
  geom_point(data = DK_FID_BOB, size=2.5, shape=21, colour="orange") +
  geom_point(data = DK_FID_NS, size=2.5, shape=21, colour="green")

# Find strips glob and change the fill color of each strip
gt <- ggplot_gtable(ggplot_build(g))
strips <- which(startsWith(gt$layout$name,'strip'))
for (s in seq_along(strips)) {gt$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- fill_colors[s]}

# Make the plot
plot(gt)
