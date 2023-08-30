#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DDAP2 = DDAP2_facetPlotAllMorphospaces
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(phytools)

# inputs
PCA <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/PCAcomps_OUTPUTS.csv") # csv of PCA data
new.tree <- read.tree("~/Notothenioids_research/repository/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/new.tree.tre") # import trimmed phylogenetic tree

depth.0TO100 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.0TO100_OUTPUTS.RData")
depth.101TO200 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.101TO200_OUTPUTS.RData")
depth.201TO300 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.201TO300_OUTPUTS.RData")
depth.301TO400 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.301TO400_OUTPUTS.RData")
depth.401TO500 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.401TO500_OUTPUTS.RData")
depth.501TO600 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.501TO600_OUTPUTS.RData")
depth.601TO700 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.601TO700_OUTPUTS.RData")
depth.701TO800 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.701TO800_OUTPUTS.RData")

depth.0TO200 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.0TO200_OUTPUTS.RData")
depth.201TO400 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.201TO400_OUTPUTS.RData")
depth.401TO600 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.401TO600_OUTPUTS.RData")
depth.601TO800 <- readRDS("~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/AntarcticPeninsula/DDAP1_pholomorpho+depthPlots/DDAP.depth.601TO800_OUTPUTS.RData")

## make facet plot for all morpho spaces by 100s
# make morpho space for 0-100 
positions <- which(PCA$X %in% depth.0TO100$Taxon)
tipsTOdrop <- c(PCA$X)
diff <- setdiff(tipsTOdrop,PCA$X[positions])
temp.tree <- drop.tip(new.tree,diff)
space1_input <- data.frame(PCA$X[positions])
space1 <- phylomorphospace(temp.tree,PCA$X[positions],label="radial",xlab="PCA1",ylab="PCA2")
title(main = "Depth: 0-100")
# make morpho space for 101-200
positions <- which(rownames(PCA$x) %in% depth.101TO200$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space2_input <- data.frame(PCA$x[positions,1:2])
space2 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 101-200")
# make morpho space for 201-300 
positions <- which(rownames(PCA$x) %in% depth.201TO300$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space3_input <- data.frame(PCA$x[positions,1:2])
space3 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main = "Depth: 201-300")
# make morpho space for 301-400
positions <- which(rownames(PCA$x) %in% depth.301TO400$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space4_input <- data.frame(PCA$x[positions,1:2])
space4 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 301-400")
# make morpho space for 401-500
positions <- which(rownames(PCA$x) %in% depth.401TO500$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space5_input <- data.frame(PCA$x[positions,1:2])
space5 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 401-500")
# make morpho space for 501-600
positions <- which(rownames(PCA$x) %in% depth.501TO600$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space6_input <- data.frame(PCA$x[positions,1:2])
space6 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 501-600")
# make morpho space for 601-700
positions <- which(rownames(PCA$x) %in% depth.601TO700$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space7_input <- data.frame(PCA$x[positions,1:2])
space7 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 601-700")
# make morpho space for 701-800
positions <- which(rownames(PCA$x) %in% depth.701TO800$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space8_input <- data.frame(PCA$x[positions,1:2])
space8 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 701-800")

# create a data frame of the plots with the phytools components
space1_df <- phytools2ggplot(space1)
space2_df <- phytools2ggplot(space2)
space3_df <- phytools2ggplot(space3)
space4_df <- phytools2ggplot(space4)
space5_df <- phytools2ggplot(space5)
space6_df <- phytools2ggplot(space6)
space7_df <- phytools2ggplot(space7)
space8_df <- phytools2ggplot(space8)

# create a list of the plots
plots1 <- list(
  ggplot(data = space1_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space1_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 0-100")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space2_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space2_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 101-200")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space3_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space3_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 201-300")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space4_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space4_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 301-400")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space5_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space5_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 401-500")
  + xlim(-0.2, 0.1) 
  + ylim(-0.1,0.1),
  ggplot(data = space6_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space6_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 501-600")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space7_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space7_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 601-700")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space8_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space8_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 701-800")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1))

# plot the list of 100s plots as a facet plot
dev.new()
grid.arrange(grobs = plots1, ncol = 2, respect = TRUE)

## make facet plot for all morpho spaces by 200s
# make morpho space for 0-200 
positions <- which(rownames(PCA$x) %in% depth.0TO200$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space9_input <- data.frame(PCA$x[positions,1:2])
space9 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main = "Depth: 0-200")
# make morpho space for 201-400
positions <- which(rownames(PCA$x) %in% depth.201TO400$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space10_input <- data.frame(PCA$x[positions,1:2])
space10 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 201-400")
# make morpho space for 401-600 
positions <- which(rownames(PCA$x) %in% depth.401TO600$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space11_input <- data.frame(PCA$x[positions,1:2])
space11 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main = "Depth: 401-600")
# make morpho space for 601-800
positions <- which(rownames(PCA$x) %in% depth.601TO800$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space12_input <- data.frame(PCA$x[positions,1:2])
space12 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main = "Depth: 601-800")

# create a dataframe of the plots with the phytools components
space9_df <- phytools2ggplot(space9)
space10_df <- phytools2ggplot(space10)
space11_df <- phytools2ggplot(space11)
space12_df <- phytools2ggplot(space12)

# create a list of the plots
plots2 <- list(
  ggplot(data = space9_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space9_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 0-200")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space10_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space10_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 201-400")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space11_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space11_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 401-600")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space12_input, aes(x = Comp1, y = Comp2)) 
  + geom_point(size = 2, color = "darkblue") 
  + geom_segment(data = space12_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") 
  + theme_classic() 
  + ggtitle("Depth: 601-800")
  + xlim(-0.2, 0.1) + ylim(-0.1,0.1))

# plot the list of 200s plots as a facet plot
dev.new()
grid.arrange(grobs = plots2, ncol = 2)

## ----------outputs----------
### PDF
pdf(file= "DDAP2.facetPhyloSpace100s+facetPhyloSpace200s_OUTPUTS.pdf")
grid.arrange(grobs = plots1, ncol = 2, respect = TRUE)
grid.arrange(grobs = plots2, ncol = 2)
dev.off() 
