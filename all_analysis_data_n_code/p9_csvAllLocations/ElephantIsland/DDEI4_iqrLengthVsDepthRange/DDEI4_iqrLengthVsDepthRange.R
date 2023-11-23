#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DDEI3 = DDEI3_iqrLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(ggplot2)

# inputs
source("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
list.of.depth.data.100s <- readRDS("~/Notothenioids-Research/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI1_pholomorpho+depthPlots/DDEI.listOfDepthData100s_OUTPUTS.RData")
list.of.depth.data.200s <- readRDS("~/Notothenioids-Research/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI1_pholomorpho+depthPlots/DDEI.listOfDepthData200s_OUTPUTS.RData")

# apply function to all the depth ranges by 100s
iqr.100s <- lapply(list.of.depth.data.100s,calc.inter.quartile)
iqr.100s <- data.frame(iqr.100s)
colnames(iqr.100s) <- list("0to100","101to200","201to300","301to400")
iqr.100s <- t(iqr.100s)
iqr.100s <- data.frame(iqr.100s)
iqr.100s <- tibble::rownames_to_column(iqr.100s, "x")
colnames(iqr.100s)[2] <- "y"

# plot 100s IQR
dev.new()
ggplot(iqr.100s,aes(x = x , y = y))+ 
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2)

# apply function to all the depth ranges by 200s
iqr.200s <- lapply(list.of.depth.data.200s,calc.inter.quartile)
iqr.200s <- data.frame(iqr.200s)
colnames(iqr.200s) <- list("0to200","201to400")
iqr.200s <- t(iqr.200s)
iqr.200s <- data.frame(iqr.200s)
iqr.200s <- tibble::rownames_to_column(iqr.200s, "x")
colnames(iqr.200s)[2] <- "y"

# plot 200s IQR
dev.new()
ggplot(iqr.200s,aes(x = x , y = y))+
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+ 
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 

## Outputs
### PDF
pdf(file= "DDEI4.plot100sIQR+plot200sIQR_OUTPUTS.pdf")
ggplot(iqr.100s,aes(x = x , y = y))+ 
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 

ggplot(iqr.200s,aes(x = x , y = y))+
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+ 
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2)  
dev.off() 