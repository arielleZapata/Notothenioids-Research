#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT8 = p8_speciesVsDepth+speciesVsLength
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(ggplot2)
library(ggridges)

# inputs
depth.data <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/CombinedAntarcticData.csv") # depth data for fish
depth.data.filtered <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p2_pholomorpho+depthPlots/depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data
fishBase.nameChange <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/fishBase_names.csv") # csv of name replacement for retrieving the fishbase data
fishB.names <- fishBase.nameChange

# species vs depth
dev.new()
ggplot(depth.data,aes(x= Depth,y=Taxon, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Species by Depth") + 
  geom_jitter(width = .05, alpha = .3)

# species vs length
depth.data <- na.omit(depth.data)
outliers <- c()
for (k in 1:nrow(depth.data)){
  if (depth.data$Taxon[k] %in% fishB.names$X){
    equalNames <- which(depth.data$Taxon[k]==fishB.names$X)  
    for (m in 1:length(equalNames)){
      if (depth.data$Length[k]>fishB.names$Size[m]){
        outliers[k] <- "Yes"
      }
      else{
        outliers[k] <- "No"
      }
    }
  }
  else{
    outliers[k] <- "No"
  }
}
results <- outliers
depth.data$Outliers <- results
no.outlier.depth.data <- subset(depth.data, Outliers == "No")
dev.new()
ggplot(no.outlier.depth.data,aes(x=Length,y=Taxon, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Length", option = "C") + 
  labs(title = "Length Distribution of Species") + 
  geom_jitter(width = .05, alpha = .3)

# just test if they really delete
lepi_squa.df <- which(no.outlier.depth.data$Taxon == "lepi_squa-trematominae")
lepi.location <- which(fishB.names$X == "lepi_squa-trematominae")
size <- fishB.names[lepi.location,]$Size
lepi_squa.df <- no.outlier.depth.data[lepi_squa.df,]
new.df <- which(lepi_squa.df$Length>size)
lepi_squa.df <- lepi_squa.df[new.df,]


# check why this does not work
dev.new()
ggplot(depth.data,aes(x=Length,y=Taxon, fill = group)) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Length", option = "C") + 
  labs(title = "Length of species") + 
  geom_jitter(width = .05, alpha = .3)

## ----------outputs----------
# PDF
pdf(file= "speciesVsDepth+speciesVsLength_OUTPUTS.pdf" )
ggplot(depth.data,aes(x= Depth,y=Taxon, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + 
  scale_fill_viridis_c(name = "Depth", option = "C") + 
  labs(title = "Depth of species") + 
  geom_jitter(width = .05, alpha = .3)
dev.off() 