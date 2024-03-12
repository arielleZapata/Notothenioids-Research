#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT1 = p1_gpa_and_phylomorpho_of_all_data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(geomorph)
library(phytools)

#set wd to folder with needed content
setwd("~/Notothenioids-Research/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data")

# make variable for file locations
tps.file <- "~/Notothenioids-Research/all_analysis_data_n_code/make_tps_files+combined_p_file/all_tps_file.tps" 

# read in data
fish <- read.tree("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/notothenioid_timetree.tre") # given phylo tree
csv.names <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/updated_names.csv",header=FALSE) # csv of name replacement for tree

# read the TPS file
data <- readland.tps(tps.file,specID = c("None", "ID", "imageID"), negNA = FALSE,readcurves = FALSE, warnmsg = FALSE)

# perform and display GPA on TPS file
data.gpa <- gpagen(data,print.progress = FALSE)
summary(data.gpa)
dev.new()
plotAllSpecimens(data.gpa$coords)
title("GPA Consensus")

# read in and display time tree
fish.phylo <- as.phylo(fish)
dev.new()
plot(fish.phylo)
title("Given Phylogenetic Tree for Notothenioids")

# rename the tree to match data
tips <- c()
for(i in 1:nrow(csv.names)) {
  if (any(csv.names[i, 2]=="")){
    tips <- c(tips,csv.names[i,1])
  }
  else{
    fish.phylo$tip.label[fish.phylo$tip.label==csv.names[i, 1]] <- csv.names[i, 2]
  }
}
new.tree <- drop.tip(fish.phylo,tips)
dev.new()
plot(new.tree)
title("Pruned Given Phylogenetic Tree for Notothenioids")

# collect new tip labels for filtering the GPA data
tip.names <- c()
for (tip in new.tree$tip.label){
  tip.names <- c(tip.names,tip)
}
tip.names <- rev(tip.names)

# replaces names in GPA
dimnames(data.gpa$coords)[[3]]<-tip.names

# perform PCA on the data
PCA <- gm.prcomp(data.gpa$coords, phy=new.tree)
summary(PCA)
dev.new()
plot(PCA$x[,1:2])
title("PCA - Comp1 vs Comp2")

# plot the data in phylomorphospace
dev.new()
phylomorphospace(new.tree,PCA$x[,1:2],label=FALSE)
title(main="Phylomorphospace Plot of All Data",font.main=3)

## ----------outputs----------
# CSV
write.csv(PCA$x[,1:2], "PCAcomps_OUTPUTS.csv")

# TRE
write.tree(new.tree, file = "new.tree.tre")

# PDF
pdf(file= "gpaCons+phylotree+pcaComp_OUTPUTS.pdf" )
plotAllSpecimens(data.gpa$coords)
title("GPA Consensus")
plot(fish.phylo)
title("Given Phylogenetic Tree for Notothenioids")
plot(new.tree)
title("Pruned Given Phylogenetic Tree for Notothenioids")
plot(PCA$x[,1:2])
title("PCA - Comp1 vs Comp2")
phylomorphospace(new.tree,PCA$x[,1:2],label=FALSE)
title(main="Phylomorphospace Plot of All Data",font.main=3)
dev.off() 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT8 = p8_speciesVsDepth+speciesVsLength
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# inputs
depth.data.filtered <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p2_pholomorpho+depthPlots/depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data
merged_data <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p7_bodyLengthVsComps/merged_data.csv") # csv of the fishbase length

dev.new()
ggplot(depth.data.filtered,aes(x= Depth,y=Species, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Species by Depth") + 
  geom_jitter(width = .05, alpha = .3)

# create a new data frame that holds all lengths less than the fishbase length
# removes outliers
smaller.than.max.df <- data.frame()
targets <- unique(depth.data.filtered$Species)
for(k in 1:length(targets)){
  target <- targets[k]
  subsample <- depth.data.filtered %>% filter(Species == target)
  locationInMergedData <- which(target %in% merged_data$X)
  
  if (length(locationInMergedData) > 0) {
    body.size <- subsample$Length
    max.length <- merged_data$SizeData[locationInMergedData]
    smaller.than.max <- subsample %>% filter(body.size < max.length)
    smaller.than.max.df <- rbind(smaller.than.max.df, smaller.than.max)
  }
}

# plot species vs filtered length
dev.new()
ggplot(smaller.than.max.df,aes(x=Length,y=Species,fill = after_stat(x))) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Length (cm)", option = "C") + 
  labs(title = "Length Distribution of Species") + 
  geom_jitter(width = .05, alpha = .3)

## ----------outputs----------
# PDF
pdf(file= "p8.speciesVsDepth+speciesVsLength_OUTPUTS.pdf" )
ggplot(depth.data.filtered,aes(x= Depth,y=Species, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Species by Depth") + 
  geom_jitter(width = .05, alpha = .3)

ggplot(smaller.than.max.df,aes(x=Length,y=Species,fill = after_stat(x))) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Length (cm)", option = "C") + 
  labs(title = "Length Distribution of Species") + 
  geom_jitter(width = .05, alpha = .3)
dev.off() 

library(viridis)
my_blue_palette <- colorRampPalette(c("lightblue", "darkblue"))

pdf(file = "p8.speciesVsDepth+speciesVsLength_OUTPUTS2.pdf")
ggplot(depth.data.filtered, aes(x = Depth, y = Species, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, palette = my_blue_palette) +
  scale_fill_gradientn(colors = my_blue_palette(100), limits = c(min(depth.data.filtered$Depth), max(depth.data.filtered$Depth)), guide = "colorbar") +
  labs(title = "Species by Depth")
dev.off()



# Assuming depth.data.filtered and fish.phylo are already loaded as shown in the original code

# Extract unique species names from the depth data
species_from_depth_data <- unique(depth.data.filtered$Species)

# Load the phylogenetic tree
fish.phylo <- read.tree("path/to/your/phylogenetic/tree/file.tre")

# Potentially prune the tree to contain only the species present in the depth data
pruned_tree <- drop.tip(fish.phylo, fish.phylo$tip.label[!(fish.phylo$tip.label %in% species_from_depth_data)])




# Example renaming, assuming you have a dataframe `name_mapping` with columns 'TreeName' and 'DepthDataName'
for(i in 1:nrow(name_mapping)) {
  pruned_tree$tip.label[pruned_tree$tip.label == name_mapping$TreeName[i]] <- name_mapping$DepthDataName[i]
}
