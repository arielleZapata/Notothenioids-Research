#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT5 = p5_iqrLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# inputs
readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p2_pholomorpho+depthPlots/list.of.depth.data.100s_OUTPUTS.rds")
readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p2_pholomorpho+depthPlots/list.of.depth.data.200s_OUTPUTS.rds")

# IQR length vs depth range
# apply function to all the depth ranges by 100s
iqr.100s <- lapply(list.of.depth.data.100s,calc.inter.quartile)
iqr.100s <- data.frame(iqr.100s)
colnames(iqr.100s) <- list("0to100","101to200","201to300","301to400","401to500","501to600","601to700","701to800")
iqr.100s <- t(iqr.100s)
iqr.100s <- data.frame(iqr.100s)
iqr.100s <- tibble::rownames_to_column(iqr.100s, "x")
colnames(iqr.100s)[2] <- "y"

# plot 100s IQR
dev.new()
ggplot(iqr.100s,aes(x = x , y = y)) 
  + geom_segment( aes(x=x, xend=x, y=0, yend=y)) 
  + geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 

# apply function to all the depth ranges by 200s
iqr.200s <- lapply(list.of.depth.data.200s,calc.inter.quartile)
iqr.200s <- data.frame(iqr.200s)
colnames(iqr.200s) <- list("0to200","201to400","401to600","601to800")
iqr.200s <- t(iqr.200s)
iqr.200s <- data.frame(iqr.200s)
iqr.200s <- tibble::rownames_to_column(iqr.200s, "x")
colnames(iqr.200s)[2] <- "y"

# plot 200s IQR
dev.new()
ggplot(iqr.200s,aes(x = x , y = y)) 
  + geom_segment( aes(x=x, xend=x, y=0, yend=y)) 
  + geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 


# Notes
## Inputs
### list.of.depth.data.100s 
### list.of.depth.data.200s


## Outputs
### PDF - plot 100s IQR -- plot 200s IQR
