#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DDEI5 = DDEI5_varianceLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# inputs
source("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
list.of.depth.data.100s <- readRDS("~/Notothenioids-Research/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI1_pholomorpho+depthPlots/DDEI.listOfDepthData100s_OUTPUTS.RData")
list.of.depth.data.200s <- readRDS("~/Notothenioids-Research/all_analysis_data_n_code/p9_csvAllLocations/ElephantIsland/DDEI1_pholomorpho+depthPlots/DDEI.listOfDepthData200s_OUTPUTS.RData")

# apply function to all the depth ranges by 100s
var.100s <- lapply(list.of.depth.data.100s,calc.variance)
var.100s <- data.frame(var.100s)
colnames(var.100s) <- list("0to100","101to200","201to300","301to400")
var.100s <- t(var.100s)
var.100s <- data.frame(var.100s)
var.100s <- tibble::rownames_to_column(var.100s, "x")
colnames(var.100s)[2] <- "y"

# plot 100s variance
dev.new()
ggplot(var.100s,aes(x = x,y = y))+ 
  geom_segment(aes(x=x, xend=x, y=0, yend=y))+ 
  geom_point(size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 

# apply function to all the depth ranges by 200s
var.200s <- lapply(list.of.depth.data.200s,calc.variance)
var.200s <- data.frame(var.200s)
colnames(var.200s) <- list("0to200","201to400")
var.200s <- t(var.200s)
var.200s <- data.frame(var.200s)
var.200s <- tibble::rownames_to_column(var.200s, "x")
colnames(var.200s)[2] <- "y"

# plot 200s variance
dev.new()
ggplot(var.200s,aes(x = x,y = y))+  
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+  
  geom_point( size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 

## Outputs
### PDF
pdf(file= "DDEI5.plot100svariance+plot200svariance_OUTPUTS.pdf")
ggplot(var.100s,aes(x = x,y = y))+ 
  geom_segment(aes(x=x, xend=x, y=0, yend=y))+ 
  geom_point(size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 

ggplot(var.200s,aes(x = x,y = y))+  
  geom_segment( aes(x=x, xend=x, y=0, yend=y))+  
  geom_point( size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 
dev.off() 