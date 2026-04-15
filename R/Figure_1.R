## FIGURE 1 ##

library('ncdf4')
library('marmap')
library('gdistance')
library('tidyverse')
library('Matrix')
library('igraph')
library('ggOceanMaps')
library('ggnewscale')
library('ggrepel')
library('ggpubr')
library('ggspatial')

setwd("../data")

# Load sample release coordinates
RELEAS <- read.table("SAMPLE_RELEASE_COORDS.txt", header=T)
POPX <- aggregate(x=RELEAS$X,by=list(RELEAS$population),FUN = mean)
colnames(POPX) <- c("POP","X")
POPY <- aggregate(x=RELEAS$Y,by=list(RELEAS$population),FUN = mean)
colnames(POPY) <- c("POP","Y")
POP_COORDS <- left_join(POPX,POPY,by='POP')
POP_COORDS <- POP_COORDS[c(2,7,6,5,1,10,9,3,8,4),]

# FIG 1A # map version with the ICES polygons
myPal <- colorRampPalette(c("red","gold","blue"))

P1 = basemap(limits = c(-12,7,42,54), bathymetry = TRUE,  land.col = "#eeeac4", grid.col = NA, grid.size = 0.5, base_size = 18, bathy.style = "poly_greys") +
  labs(x = "Longitude", y = "Latitude") +
  ggnewscale::new_scale_fill() +
  annotation_spatial(ices_areas[36:37,], mapping = aes(), fill="purple", col="purple", alpha = 0.1) + 
  annotation_spatial(ices_areas[39:40,], mapping = aes(), fill="orange", col="orange", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(7:8,43:48),], mapping = aes(), fill="green", col="green", alpha = 0.1) + 
  annotation_spatial(ices_areas[c(27,41:42),], mapping = aes(), fill="blue", col="blue", alpha = 0.1) + 
  theme(legend.position = "none")

P1 = reorder_layers(P1) + annotation_scale(location = "br") +
  geom_sf_label(data = suppressWarnings(sf::st_centroid(sf::st_make_valid(ices_areas[36:37,]))), col="purple", aes(label = Area_Full), fill=NA, label.size = NA) +
  geom_sf_label(data = suppressWarnings(sf::st_centroid(sf::st_make_valid(ices_areas[39:40,]))), col="orange", aes(label = Area_Full), fill=NA, label.size = NA) +
  geom_sf_label(data = suppressWarnings(sf::st_centroid(sf::st_make_valid(ices_areas[c(7:8,43:48),]))), col="green", aes(label = Area_Full), fill=NA, label.size = NA) +
  geom_sf_label(data = suppressWarnings(sf::st_centroid(sf::st_make_valid(ices_areas[c(27,41:42),]))), col="blue", aes(label = Area_Full), fill=NA, label.size = NA) +
  geom_spatial_point(POP_COORDS, mapping=aes(x=X,y=Y), colour=myPal(10), cex=4, crs = 4326) +
  stat_spatial_identity(POP_COORDS, mapping=aes(x=X,y=Y, label = POP), colour=myPal(10), geom = "label_repel", cex=5) + ggtitle("Locations of tagging surveys")


# Load sample migration phenotypes
PHENO = read.table("SAMPLE_MIGRATION_PHENOTYPE.txt", header=T)

# FIG 1B # 
# Proportion of long-distance migrants per location  
MIG_PROP <- PHENO %>% count(population,mvt_strategies) %>% na.exclude() %>% group_by(population) %>% mutate(freq = 100*(n/sum(n))) %>% filter(mvt_strategies=="migrant")
MIG_PROP <- MIG_PROP[c(2,7,6,5,1,10,9,3,8,4),]
MIG_PROP$population <- factor(MIG_PROP$population, levels = MIG_PROP$population)
myPal <- colorRampPalette(c("red","gold","blue"))
P2 = ggplot(MIG_PROP, aes(x = population, y = freq)) + geom_bar(stat="identity",fill=myPal(10)) + labs(x = "", y = "Long-distance migrants (%)") + theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=20))

# FIG 1C # 
# Proportion of spawning fidelity per location  
FID_PROP <- PHENO %>% count(population,spawning_behavior) %>% na.exclude() %>% group_by(population) %>% mutate(freq = 100*(n/sum(n))) %>% filter(spawning_behavior=="Fidelity")
FID_PROP <- FID_PROP[c(2,7,6,5,1,10,9,3,8,4),]
FID_PROP$population <- factor(FID_PROP$population, levels = FID_PROP$population)
myPal <- colorRampPalette(c("red","gold","blue"))
P3 = ggplot(FID_PROP, aes(x = population, y = freq)) + geom_bar(stat="identity",fill=myPal(10)) + labs(x = "", y = "Spawning site fidelity (%)") + theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=20))

# FIG 1D # 
# Proportion of NS/BOB fidelity per location  
NSBOB_PROP <- PHENO %>% count(population,spawning_ground_strategy) %>% na.exclude() %>% group_by(population) %>% mutate(freq = 100*(n/sum(n))) %>% group_by(spawning_ground_strategy)
NSBOB_PROP <- NSBOB_PROP[c(3,10,9,8,1,2,13,12,4,5,11,6,7),]
P4 = ggplot(NSBOB_PROP, aes(fill = spawning_ground_strategy, x = population, y = freq)) + geom_bar(stat="identity", show.legend = F) + labs(x = "", y = "NS/BOB spawning site fidelity") + theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=20)) + scale_x_discrete(limits=rev(POPS)) +
  scale_fill_manual(values = c("orange","green"))

P5 = ggarrange(P2, P3, P4, labels = c("B","C","D"), ncol = 3, nrow = 1, common.legend = F)


# Load sample summer and winter coordinates (07/09/22)
COORDS <- read.table("SAMPLE_SUMMER_WINTER_COORDS.txt", header=T)
# Sort individuals by sample name 
COORDS <- arrange(COORDS, sample_filename)

# FIG 1E # Map the SUMMER1 Coordinates
myPal <- colorRampPalette(c("blue","gold","red"))
POPS = c("DK","PB","CH","SM","SQ","AD","LT","OL","NO","CB")

IND_COORDS = NULL
IND_COL = NULL
for (p in 1:length(POPS))
  {
  X <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Summer1_X,Summer1_Y))$Summer1_X
  Y <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Summer1_X,Summer1_Y))$Summer1_Y
  POP_P <- cbind(X,Y)
  COL <- rep(myPal(10)[p],nrow(POP_P))
  IND_COORDS <- rbind(IND_COORDS,POP_P)
  IND_COL <- c(IND_COL,COL)
  }
IND_COORDS = data.frame(IND_COORDS)

P6 = basemap(limits = c(-12,7,42,54), bathymetry = TRUE,  land.col = "#eeeac4", grid.col = NA, grid.size = 0.5, base_size = 18, bathy.style = "poly_greys") +
  labs(x = "Longitude", y = "Latitude") +
  ggnewscale::new_scale_fill() +
  annotation_spatial(ices_areas[36:37,], mapping = aes(), fill="purple", col="purple", alpha = 0) + 
  annotation_spatial(ices_areas[39:40,], mapping = aes(), fill="orange", col="orange", alpha = 0) + 
  annotation_spatial(ices_areas[c(7:8,43:48),], mapping = aes(), fill="green", col="green", alpha = 0) + 
  annotation_spatial(ices_areas[c(27,41:42),], mapping = aes(), fill="blue", col="blue", alpha = 0) + 
  theme(legend.position = "none")
P6 = reorder_layers(P6) + annotation_scale(location = "br") +
  geom_spatial_point(IND_COORDS, mapping=aes(x=X,y=Y), colour=IND_COL, cex=2, crs = 4326) + ggtitle("Individuals positions - First Summer")


# FIG 1F # Map the WINTER1 Coordinates
myPal <- colorRampPalette(c("blue","gold","red"))
POPS = c("DK","PB","CH","SM","SQ","AD","LT","OL","NO","CB")

IND_COORDS = NULL
IND_COL = NULL
for (p in 1:length(POPS))
  {
  X <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_X
  Y <- (COORDS %>% filter(population==POPS[p]) %>% drop_na(Winter1_X,Winter1_Y))$Winter1_Y
  POP_P <- cbind(X,Y)
  COL <- rep(myPal(10)[p],nrow(POP_P))
  IND_COORDS <- rbind(IND_COORDS,POP_P)
  IND_COL <- c(IND_COL,COL)
  }
IND_COORDS = data.frame(IND_COORDS)

P7 = basemap(limits = c(-12,7,42,54), bathymetry = TRUE,  land.col = "#eeeac4", grid.col = NA, grid.size = 0.5, base_size = 18, bathy.style = "poly_greys") +
  labs(x = "Longitude", y = "Latitude") +
  ggnewscale::new_scale_fill() +
  annotation_spatial(ices_areas[36:37,], mapping = aes(), fill="purple", col="purple", alpha = 0) + 
  annotation_spatial(ices_areas[39:40,], mapping = aes(), fill="orange", col="orange", alpha = 0) + 
  annotation_spatial(ices_areas[c(7:8,43:48),], mapping = aes(), fill="green", col="green", alpha = 0) + 
  annotation_spatial(ices_areas[c(27,41:42),], mapping = aes(), fill="blue", col="blue", alpha = 0) + 
  theme(legend.position = "none")
P7 = reorder_layers(P7) + annotation_scale(location = "br") +
  geom_spatial_point(IND_COORDS, mapping=aes(x=X,y=Y), colour=IND_COL, cex=2, crs = 4326) + ggtitle("Individuals positions - First Winter")


# PLOT WHOLE FIG 1 #
P8 = ggarrange(P1, P5, P6, P7, labels = c("A","","E","F"), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
P8
