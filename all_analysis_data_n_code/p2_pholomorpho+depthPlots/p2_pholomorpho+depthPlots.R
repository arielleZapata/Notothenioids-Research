#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT2 = p2_phylomorpho_of_all_data+depthDataSubsets+convexHullGraphs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(dplyr)
library(phytools)
library(geiger)
library(splancs)

# inputs
source("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
depth.data <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/CombinedAntarcticData.csv") # depth data for fish
depthNameChange <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/FIXEDdepthNames.csv") # csv of name replacement for depth data
PCA <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/PCAcomps_OUTPUTS.csv") # csv of PCA data
new.tree <- read.tree("~/Notothenioids_research/repository/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/new.tree.tre") # import trimmed phylogenetic tree

# location site input
depth.data.Antarctic.Peninsula <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.Antarctic.Peninsula.csv") # depth data for Antarctic Peninsula fish
depth.data.Elephant.Island <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.Elephant.Island.csv") # depth data for Antarctic Peninsula fish
depth.data.South.Orkney.Is <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.South.Orkney.Is.csv") # depth data for Antarctic Peninsula fish
depth.data.South.Shetland.Is <-  read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.South.Shetland.Is.csv") # depth data for Antarctic Peninsula fish

# change depth data names to names in tree
for (i in 1:length(depth.data$Taxon)){
  for (j in 1:length(depthNameChange$Old)){
    if (depth.data$Taxon[i] == depthNameChange$Old[j]){
      depth.data$Taxon[i]<-depthNameChange$New[j]
    }
  }
}

# filter ~ missing names from depth data
depth.data <- depth.data %>% filter(!Taxon=="~")

# subset data using 0-100,101-200,201-300,301-400,401-500,501-600,601-700,701-800
depth.0TO100 <- subset(depth.data,depth.data$Depth>0 & depth.data$Depth<=100)
depth.101TO200 <- subset(depth.data,depth.data$Depth>101 & depth.data$Depth<=200)
depth.201TO300 <- subset(depth.data,depth.data$Depth>201 & depth.data$Depth<=300)
depth.301TO400 <- subset(depth.data,depth.data$Depth>301 & depth.data$Depth<=400)
depth.401TO500 <- subset(depth.data,depth.data$Depth>401 & depth.data$Depth<=500)
depth.501TO600 <- subset(depth.data,depth.data$Depth>501 & depth.data$Depth<=600)
depth.601TO700 <- subset(depth.data,depth.data$Depth>601 & depth.data$Depth<=700)
depth.701TO800 <- subset(depth.data,depth.data$Depth>701 & depth.data$Depth<=800)

# subset data using 0-200,201-400,401-600,601-800
depth.0TO200 <- subset(depth.data,depth.data$Depth>0 & depth.data$Depth<=200)
depth.201TO400 <- subset(depth.data,depth.data$Depth>201 & depth.data$Depth<=400)
depth.401TO600 <- subset(depth.data,depth.data$Depth>401 & depth.data$Depth<=600)
depth.601TO800 <- subset(depth.data,depth.data$Depth>601 & depth.data$Depth<=800)

# create new df based on the species from phylomorphospace
depth.data.filtered <- rbind(depth.0TO100,depth.101TO200,depth.201TO300,depth.301TO400,depth.401TO500,depth.501TO600,depth.601TO700,depth.701TO800)
list.of.depth.data.100s <- list(depth.0TO100,depth.101TO200,depth.201TO300,depth.301TO400,depth.401TO500,depth.501TO600,depth.601TO700,depth.701TO800)
list.of.depth.data.200s <- list(depth.0TO200,depth.201TO400,depth.401TO600,depth.601TO800)

# run function on 0-800 by 100s
a1 <- calcConvex.phylomorpho.chull(depth.0TO100,"Depth: 0-100",new.tree)
a2 <- calcConvex.phylomorpho.chull(depth.101TO200,"Depth: 101-200",new.tree)
a3 <- calcConvex.phylomorpho.chull(depth.201TO300,"Depth: 201-300",new.tree)
a4 <- calcConvex.phylomorpho.chull(depth.301TO400,"Depth: 301-400",new.tree)
a5 <- calcConvex.phylomorpho.chull(depth.401TO500,"Depth: 401-500",new.tree)
a6 <- calcConvex.phylomorpho.chull(depth.501TO600,"Depth: 501-600",new.tree)
a7 <- calcConvex.phylomorpho.chull(depth.601TO700,"Depth: 601-700",new.tree)
a8 <- calcConvex.phylomorpho.chull(depth.701TO800,"Depth: 701-800",new.tree)

# calculate and plot the areas for the 100s
area.val <- c(a1[[2]],a2[[2]],a3[[2]],a4[[2]],a5[[2]],a6[[2]],a7[[2]],a8[[2]])
dev.new()
plot(area.val,main="Depth Range vs CHull Area by the 100s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")

# run function on 0-800 by 200s
b1 <- calcConvex.phylomorpho.chull(depth.0TO200,"Depth: 0-200",new.tree)
b2 <- calcConvex.phylomorpho.chull(depth.201TO400,"Depth: 201-400",new.tree)
b3 <- calcConvex.phylomorpho.chull(depth.401TO600,"Depth: 401-600",new.tree)
b4 <- calcConvex.phylomorpho.chull(depth.601TO800,"Depth: 601-800",new.tree)

# calculate and plot the areas for the 200s
area.val <- c(b1[[2]],b2[[2]],b3[[2]],b4[[2]])
dev.new()
plot(area.val,main="Depth Range vs CHull Area by the 200s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")

## ----------outputs----------
### CSV - depth.data.filtered
write.csv(depth.data.filtered, "depth.data.filtered_OUTPUTS.csv")

### RData
saveRDS(object = list.of.depth.data.100s, file = "listOfDepthData100s_OUTPUTS.RData")
saveRDS(object = list.of.depth.data.200s, file = "listOfDepthData200s_OUTPUTS.RData")

saveRDS(object = depth.0TO100, file = "depth.0TO100_OUTPUTS.RData")
saveRDS(object = depth.101TO200, file = "depth.101TO200_OUTPUTS.RData")
saveRDS(object = depth.201TO300, file = "depth.201TO300_OUTPUTS.RData")
saveRDS(object = depth.301TO400, file = "depth.301TO400_OUTPUTS.RData")
saveRDS(object = depth.401TO500, file = "depth.401TO500_OUTPUTS.RData")
saveRDS(object = depth.501TO600, file = "depth.501TO600_OUTPUTS.RData")
saveRDS(object = depth.601TO700, file = "depth.601TO700_OUTPUTS.RData")
saveRDS(object = depth.701TO800, file = "depth.701TO800_OUTPUTS.RData")

saveRDS(object = depth.0TO200, file = "depth.0TO200_OUTPUTS.RData")
saveRDS(object = depth.201TO400, file = "depth.201TO400_OUTPUTS.RData")
saveRDS(object = depth.401TO600, file = "depth.401TO600_OUTPUTS.RData")
saveRDS(object = depth.601TO800, file = "depth.601TO800_OUTPUTS.RData")


### PDF
pdf(file= "cHullAreas100s+cHullAreas100s_OUTPUTS.pdf")
plot(area.val,main="Depth Range vs CHull Area by the 100s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")
plot(area.val,main="Depth Range vs CHull Area by the 200s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")
dev.off() 

