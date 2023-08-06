#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT4 = p4_depthRangeVsLengthPlot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(ggplot2)
library(ggdist)
library(tidyquant)

# inputs
depth.data.filtered <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p2_pholomorpho+depthPlots/depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data

# ADD FILTERED DATA FOR EACH LOCATION

# create a vector to store the group identities
groups <- character(length(depth.data.filtered$Depth))

# assign group identities based on data values
groups[depth.data.filtered$Depth >= 0 & depth.data.filtered$Depth <= 100] <- "0to100"
groups[depth.data.filtered$Depth >= 101 & depth.data.filtered$Depth <= 200] <- "101to200"
groups[depth.data.filtered$Depth >= 201 & depth.data.filtered$Depth <= 300] <- "201to300"
groups[depth.data.filtered$Depth >= 301 & depth.data.filtered$Depth <= 400] <- "301to400"
groups[depth.data.filtered$Depth >= 401 & depth.data.filtered$Depth <= 500] <- "401to500"
groups[depth.data.filtered$Depth >= 501 & depth.data.filtered$Depth <= 600] <- "501to600"
groups[depth.data.filtered$Depth >= 601 & depth.data.filtered$Depth <= 700] <- "601to700"
groups[depth.data.filtered$Depth >= 701 & depth.data.filtered$Depth <= 800] <- "701to800"

# add the group column to the data frame
depth.data.filtered$group <- groups

# raincloud plot
dev.new()
ggplot(depth.data.filtered, aes(y = Length, x = group, fill = group)) +  
  stat_halfeye(adjust = 0.5,justification = -0.2,.width = 0,point_colour = NA) +  
  geom_boxplot(width = 0.12,outlier.color = NA,alpha = 0.5) +  
  geom_jitter(width = .05, alpha = .3) +  
  scale_fill_tq() +  
  theme_tq() +
  labs(title = "Depth Range Vs Length RainCloud Plot",x = "Depth",y = "Lengths",fill = "Depth") +
  coord_flip()

## ----------outputs----------
### PDF
pdf(file= "depthRangeVsLengthPlot_OUTPUTS.pdf")
ggplot(depth.data.filtered, aes(y = Length, x = group, fill = group)) +  
  stat_halfeye(adjust = 0.5,justification = -0.2,.width = 0,point_colour = NA) +  
  geom_boxplot(width = 0.12,outlier.color = NA,alpha = 0.5) +  
  geom_jitter(width = .05, alpha = .3) +  
  scale_fill_tq() +  
  theme_tq() +
  labs(title = "Depth Range Vs Length RainCloud Plot",x = "Depth",y = "Lengths",fill = "Depth") +
  coord_flip()
dev.off() 