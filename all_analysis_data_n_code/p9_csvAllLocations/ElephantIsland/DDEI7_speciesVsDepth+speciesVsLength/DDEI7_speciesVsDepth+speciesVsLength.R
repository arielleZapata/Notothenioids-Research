#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DDEI7 = DDEI7_speciesVsDepth+speciesVsLength
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(ggplot2)
library(ggridges)
library(dplyr)

# inputs
depth.data <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/depth.data.Elephant.Island.csv") # depth data for fish
depth.data.filtered <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI1_pholomorpho+depthPlots/DDEI.depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data
merged_data.revLog <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI6_bodyLengthVsComps/merged_data.csv") # csv of the fishbase length

dev.new()
ggplot(depth.data,aes(x= Depth,y=Taxon, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Species by Depth") + 
  geom_jitter(width = .05, alpha = .3)

# create a new data frame that holds all lengths less than the fishbase length
# removes outliers
#smaller.than.max.df <- data.frame()
#targets <- unique(depth.data.filtered$Taxon)
#for(k in 1:length(targets)){
#  target <- targets[k]
#  subsample <- depth.data.filtered %>% filter(Taxon == target)
#  locationInMergedData <- which(target %in% merged_data$X)
  
#  if (length(locationInMergedData) > 0) {
#    body.size <- subsample$Length
#    max.length <- merged_data$SizeData[locationInMergedData]
#    smaller.than.max <- subsample %>% filter(body.size < max.length)
#    smaller.than.max.df <- rbind(smaller.than.max.df, smaller.than.max)
#  }
#}

# plot species vs filtered length
#dev.new()
#ggplot(smaller.than.max.df,aes(x=Length,y=Taxon,fill = after_stat(x))) + 
#  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
#  scale_fill_viridis_c(name = "Length (cm)", option = "C") + 
#  labs(title = "Length Distribution of Species") + 
#  geom_jitter(width = .05, alpha = .3)

## ----------outputs----------
# PDF
pdf(file= "DDEI7.speciesVsDepth+speciesVsLength_OUTPUTS.pdf" )
ggplot(depth.data,aes(x= Depth,y=Taxon, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Species by Depth") + 
  geom_jitter(width = .05, alpha = .3)

#ggplot(smaller.than.max.df,aes(x=Length,y=Taxon,fill = after_stat(x))) + 
#  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
#  scale_fill_viridis_c(name = "Length (cm)", option = "C") + 
#  labs(title = "Length Distribution of Species") + 
#  geom_jitter(width = .05, alpha = .3)
dev.off() 