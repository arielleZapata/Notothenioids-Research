#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT3 = p3_facetPlotAllMorphospaces
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(ggplot2)
library(phytools)
library(gridExtra)

# inputs
source("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
setwd("~/Notothenioids-Research/all_analysis_data_n_code/p2_pholomorpho+depthPlots") # location of all the R lists data
PCA <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/pca.csv") # csv of PCA data
new.tree <- read.tree("~/Notothenioids-Research/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/new.tree.tre") # import trimmed phylogenetic tree

depth.0TO100 <- readRDS("depth.0TO100_OUTPUTS.RData")
depth.101TO200 <- readRDS("depth.101TO200_OUTPUTS.RData")
depth.201TO300 <- readRDS("depth.201TO300_OUTPUTS.RData")
depth.301TO400 <- readRDS("depth.301TO400_OUTPUTS.RData")
depth.401TO500 <- readRDS("depth.401TO500_OUTPUTS.RData")
depth.501TO600 <- readRDS("depth.501TO600_OUTPUTS.RData")
depth.601TO700 <- readRDS("depth.601TO700_OUTPUTS.RData")
depth.701TO800 <- readRDS("depth.701TO800_OUTPUTS.RData")
depth.greaterThan801 <- readRDS("depth.greaterThan801_OUTPUTS.RData")

depth.0TO200 <- readRDS("depth.0TO200_OUTPUTS.RData")
depth.201TO400 <- readRDS("depth.201TO400_OUTPUTS.RData")
depth.401TO600 <- readRDS("depth.401TO600_OUTPUTS.RData")
depth.601TO800 <- readRDS("depth.601TO800_OUTPUTS.RData")

## make facet plot for all morpho spaces by 100s
space1_list <- make.morpho.depth(PCA,depth.0TO100,new.tree,title="Depth: 0-100")
space2_list <-make.morpho.depth(PCA,depth.101TO200,new.tree,title="Depth: 101-200")
space3_list <-make.morpho.depth(PCA,depth.201TO300,new.tree,title="Depth: 201-300")
space4_list <-make.morpho.depth(PCA,depth.301TO400,new.tree,title="Depth: 301-400")
space5_list <-make.morpho.depth(PCA,depth.401TO500,new.tree,title="Depth: 401-500")
space6_list <-make.morpho.depth(PCA,depth.501TO600,new.tree,title="Depth: 501-600")
space7_list <-make.morpho.depth(PCA,depth.601TO700,new.tree,title="Depth: 601-700")
space8_list <-make.morpho.depth(PCA,depth.701TO800,new.tree,title="Depth: 701-800")
space9_list <-make.morpho.depth(PCA,depth.greaterThan801,new.tree,title="Depth: 800+",dfName="Species")


# create a data frame of the plots with the phytools components
space1_df <- phytools2ggplot(space1_list[[2]])
space2_df <- phytools2ggplot(space2_list[[2]])
space3_df <- phytools2ggplot(space3_list[[2]])
space4_df <- phytools2ggplot(space4_list[[2]])
space5_df <- phytools2ggplot(space5_list[[2]])
space6_df <- phytools2ggplot(space6_list[[2]])
space7_df <- phytools2ggplot(space7_list[[2]])
space8_df <- phytools2ggplot(space8_list[[2]])
space9_df <- phytools2ggplot(space9_list[[2]])

# create a list of the plots
plots1 <- list(
  ggplot(data = space1_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue")+ 
    geom_segment(data = space1_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 0-100")+ 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space2_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space2_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 101-200") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space3_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space3_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 201-300") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space4_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space4_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 301-400") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space5_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space5_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 401-500") + 
    xlim(-0.2, 0.1) + 
    ylim(-0.1,0.1),
  ggplot(data = space6_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space6_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 501-600") + 
    xlim(-0.2, 0.1) + 
    ylim(-0.1,0.1),
  ggplot(data = space7_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space7_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 601-700") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space8_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space8_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 701-800") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space9_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space9_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 800+") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1))

# plot the list of 100s plots as a facet plot
dev.new()
grid.arrange(grobs = plots1, ncol = 2, respect = TRUE)

## make facet plot for all morpho spaces by 200s
space9_list <- make.morpho.depth(PCA,depth.0TO200,new.tree,title="Depth: 0-200")
space10_list <-make.morpho.depth(PCA,depth.201TO400,new.tree,title="Depth: 201-400")
space11_list <-make.morpho.depth(PCA,depth.401TO600,new.tree,title="Depth: 401-600")
space12_list <-make.morpho.depth(PCA,depth.601TO800,new.tree,title="Depth: 601-800")
space13_list <-make.morpho.depth(PCA,depth.greaterThan801,new.tree,title="Depth: 800+",dfName="Species")

# create a dataframe of the plots with the phytools components
space9_df <- phytools2ggplot(space9_list[[2]])
space10_df <- phytools2ggplot(space10_list[[2]])
space11_df <- phytools2ggplot(space11_list[[2]])
space12_df <- phytools2ggplot(space12_list[[2]])
space13_df <- phytools2ggplot(space13_list[[2]])

# create a list of the plots
plots2 <- list(
  ggplot(data = space9_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space9_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") +
    theme_classic() + 
    ggtitle("Depth: 0-200") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space10_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space10_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 201-400") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space11_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space11_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 401-600") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space12_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space12_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 601-800") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  ggplot(data = space13_list[[1]], aes(x = X1, y = X2)) + 
    geom_point(size = 2, color = "darkblue") + 
    geom_segment(data = space13_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop), linewidth = 0.5, color = "darkblue") + 
    theme_classic() + 
    ggtitle("Depth: 800+") + 
    xlim(-0.2, 0.1) + ylim(-0.1,0.1)
  )

# plot the list of 200s plots as a facet plot
dev.new()
grid.arrange(grobs = plots2, ncol = 2)

## ----------outputs----------
### PDF
setwd("~/Notothenioids-Research/all_analysis_data_n_code/p3_facetPlotAllMorphospaces")
pdf(file= "facetPhyloSpace100s+facetPhyloSpace200s_OUTPUTS.pdf")
grid.arrange(grobs = plots1, ncol = 2, respect = TRUE)
grid.arrange(grobs = plots2, ncol = 2)
dev.off() 