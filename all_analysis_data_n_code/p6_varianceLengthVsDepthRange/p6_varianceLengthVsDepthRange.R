#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT6 = p6_varianceLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# inputs
## list.of.depth.data.100s 
## list.of.depth.data.200s

# Variance length vs depth range = regular plot (lollipop)
# Calculate coefficient of variance ( sd/mean ) same as IQR graph 

# apply function to all the depth ranges by 100s
var.100s <- lapply(list.of.depth.data.100s,calc.variance)
var.100s <- data.frame(var.100s)
colnames(var.100s) <- list("0to100","101to200","201to300","301to400","401to500","501to600","601to700","701to800")
var.100s <- t(var.100s)
var.100s <- data.frame(var.100s)
var.100s <- tibble::rownames_to_column(var.100s, "x")
colnames(var.100s)[2] <- "y"
# plot 100s variance
dev.new()
ggplot(var.100s,aes(x = x,y = y)) 
  + geom_segment(aes(x=x, xend=x, y=0, yend=y)) 
  + geom_point(size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 

# apply function to all the depth ranges by 200s
var.200s <- lapply(list.of.depth.data.200s,calc.variance)
var.200s <- data.frame(var.200s)
colnames(var.200s) <- list("0to200","201to400","401to600","601to800")
var.200s <- t(var.200s)
var.200s <- data.frame(var.200s)
var.200s <- tibble::rownames_to_column(var.200s, "x")
colnames(var.200s)[2] <- "y"
# plot 200s variance
dev.new()
ggplot(var.200s,aes(x = x,y = y)) 
  +  geom_segment( aes(x=x, xend=x, y=0, yend=y)) 
  +  geom_point( size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 


# Notes
## Inputs
### list.of.depth.data.100s 
### list.of.depth.data.200s


## Outputs
### PDF - plot 100s variance -- plot 200s variance