# load in all libraries needed
library(geomorph)
library(phytools)
library(dplyr)
library(magrittr)
library(ggplot2)
library(splancs)
library(gridExtra)
library(ggdist)
library(tidyverse)
library(tidyquant)
library(ggthemes)
library(rfishbase)
library(ggridges)
library(viridis) 
library(strip)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT1 = p1_gpa_and_phylomorpho_of_all_data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# make variable for file locations
tps.file <- "~/Notothenioids_research/repository/all_analysis_data_n_code/make_tps_files+combined_p_file/all_tps_file.tps" 

# read in data
fish <- read.tree("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/notothenioid_timetree.tre") # given phylo tree
csv.names <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/updated_names.csv",header=FALSE) # csv of name replacement for tree
depth.data <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/CombinedAntarcticData.csv") # depth data for fish
depthNameChange <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/FIXEDdepthNames.csv") # csv of name replacement for depth data
fishBase.nameChange <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/fishBase_names.csv") # csv of name replacement for retrieving the fishbase data

# read the TPS file
data <- readland.tps(tps.file,specID = c("None", "ID", "imageID"), negNA = FALSE,readcurves = FALSE, warnmsg = FALSE)
print("TPS FILE COMPLETE")

# perform and display GPA on TPS file
data.gpa <- gpagen(data,print.progress = FALSE)
print("GPA ANALYSIS COMPLETE")
summary(data.gpa)
plotAllSpecimens(data.gpa$coords)

# read in and display time tree
fish.phylo <- as.phylo(fish)
print("ORIGINAL PHYLOGENETIC TREE LOADED")
plot(fish.phylo)

# rename the tree to match data
tips <- c()
for(i in 1:nrow(csv.names)) {
  if (any(csv.names[i, 2]=="")){
    tips <- c(tips,csv.names[i,1])
  }
  else{
    fish.phylo$tip.label[fish.phylo$tip.label==csv.names[i, 1]] <- csv.names[i, 2]
  }
}
new.tree <- drop.tip(fish.phylo,tips)
print("NEW PHYLOGENETIC TREE COMPLETE")
plot(new.tree)

# collect new tip labels for filtering the GPA data
tip.names <- c()
for (tip in new.tree$tip.label){
  tip.names <- c(tip.names,tip)
}
tip.names <- rev(tip.names)

# replaces names in GPA
dimnames(data.gpa$coords)[[3]]<-tip.names

# perform PCA on the data
PCA <- gm.prcomp(data.gpa$coords, phy=new.tree)
summary(PCA)
plot(PCA, main = "PCA")
plot(PCA$x[,1:2])

# plot the data in phylomorphospace
dev.new(width=5, height=2, unit="in")
phylomorphospace(new.tree,PCA$x[,1:2],label=FALSE)
title(main="Phylomorphospace Plot of All Data",font.main=3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT2 = p2_phylomorpho_of_all_data+depthDataSubsets+convexHullGraphs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (i in 1:length(depth.data$Taxon)){
  for (j in 1:length(depthNameChange$Old)){
    if (depth.data$Taxon[i] == depthNameChange$Old[j]){
      depth.data$Taxon[i]<-depthNameChange$New[j]
   }
  }
}
print("CHECK FOR PROPER NAME CHANGE")
print(unique(depth.data$Taxon))

# filter ~ missing names from depth data
depth.data<- depth.data %>% filter(!Taxon=="~")
print("CHECK TO REMOVE `~`")
print(unique(depth.data$Taxon))

# subset data using 0-100,101-200,201-300,301-400,401-500,501-600,601-700,701-800
depth.0TO100<-subset(depth.data,depth.data$Depth>0 & depth.data$Depth<=100)
depth.101TO200<-subset(depth.data,depth.data$Depth>101 & depth.data$Depth<=200)
depth.201TO300<-subset(depth.data,depth.data$Depth>201 & depth.data$Depth<=300)
depth.301TO400<-subset(depth.data,depth.data$Depth>301 & depth.data$Depth<=400)
depth.401TO500<-subset(depth.data,depth.data$Depth>401 & depth.data$Depth<=500)
depth.501TO600<-subset(depth.data,depth.data$Depth>501 & depth.data$Depth<=600)
depth.601TO700<-subset(depth.data,depth.data$Depth>601 & depth.data$Depth<=700)
depth.701TO800<-subset(depth.data,depth.data$Depth>701 & depth.data$Depth<=800)

# subset data using 0-200,201-400,401-600,601-800
depth.0TO200<-subset(depth.data,depth.data$Depth>0 & depth.data$Depth<=200)
depth.201TO400<-subset(depth.data,depth.data$Depth>201 & depth.data$Depth<=400)
depth.401TO600<-subset(depth.data,depth.data$Depth>401 & depth.data$Depth<=600)
depth.601TO800<-subset(depth.data,depth.data$Depth>601 & depth.data$Depth<=800)

# make morpho spaces + calculate convex hulls function
calcConvex.phylomorpho.chull <- function(df,title){
  positions <- which(rownames(PCA$x) %in% df$Taxon)
  tipsTOdrop <- c(rownames(PCA$x))
  diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
  temp.tree <- drop.tip(new.tree,diff)
  dev.new(width=5, height=2, unit="in")
  space<-phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
  title(main=title)
  pc <- as_tibble(PCA$x[positions,1:2])
  chull<-chull(pc)
  convex_hull <- pc %>% slice(chull)
  area <- as.matrix(convex_hull) %>% areapl
  output<-list(space,area)
  return(output)
}

# run function on 0-800 by 100s
a1 <- calcConvex.phylomorpho.chull(depth.0TO100,"Depth: 0-100")
a2 <- calcConvex.phylomorpho.chull(depth.101TO200,"Depth: 101-200")
a3 <- calcConvex.phylomorpho.chull(depth.201TO300,"Depth: 201-300")
a4 <- calcConvex.phylomorpho.chull(depth.301TO400,"Depth: 301-400")
a5 <- calcConvex.phylomorpho.chull(depth.401TO500,"Depth: 401-500")
a6 <- calcConvex.phylomorpho.chull(depth.501TO600,"Depth: 501-600")
a7 <- calcConvex.phylomorpho.chull(depth.601TO700,"Depth: 601-700")
a8 <- calcConvex.phylomorpho.chull(depth.701TO800,"Depth: 701-800")

# calculate and plot the areas for the 100s
area.val <- c(a1[[2]],a2[[2]],a3[[2]],a4[[2]],a5[[2]],a6[[2]],a7[[2]],a8[[2]])
dev.new(width=5, height=2, unit="in")
plot(area.val,main="Depth Range vs CHull Area by the 100s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")

# run function on 0-800 by 200s
b1 <- calcConvex.phylomorpho.chull(depth.0TO200,"Depth: 0-200")
b2 <- calcConvex.phylomorpho.chull(depth.201TO400,"Depth: 201-400")
b3 <- calcConvex.phylomorpho.chull(depth.401TO600,"Depth: 401-600")
b4 <- calcConvex.phylomorpho.chull(depth.601TO800,"Depth: 601-800")

# calculate and plot the areas for the 200s
area.val <- c(b1[[2]],b2[[2]],b3[[2]],b4[[2]])
dev.new(width=5, height=2, unit="in")
plot(area.val,main="Depth Range vs CHull Area by the 200s",xlab="Depth Range (ft)",ylab="CHull Area", type = "l")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT3 = p3_facetPlotAllMorphospaces
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# make facet plot for all morpho spaces by 100s
#-------------------------------------------------------------------------------
# phylomorpho to dataframe function
phytools2ggplot<-function(phylomorphospace){
  output <-data.frame(
    xstart= phylomorphospace$xx[phylomorphospace$edge[,1]],
    ystart=phylomorphospace$yy[phylomorphospace$edge[,1]],
    xstop=phylomorphospace$xx[phylomorphospace$edge[,2]],
    ystop=phylomorphospace$yy[phylomorphospace$edge[,2]],
    nodestart=phylomorphospace$edge[,1],
    nodestop=phylomorphospace$edge[,2])
  return(output)
}

# make morpho space for 0-100 
positions <- which(rownames(PCA$x) %in% depth.0TO100$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space1_input<-data.frame(PCA$x[positions,1:2])
space1 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
# ARE THESE NEEDED???
title(main="Depth: 0-100")

# make morpho space for 101-200
positions <- which(rownames(PCA$x) %in% depth.101TO200$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space2_input<-data.frame(PCA$x[positions,1:2])
space2 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 101-200")

# make morpho space for 201-300 
positions <- which(rownames(PCA$x) %in% depth.201TO300$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space3_input<-data.frame(PCA$x[positions,1:2])
space3 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main="Depth: 201-300")

# make morpho space for 301-400
positions <- which(rownames(PCA$x) %in% depth.301TO400$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space4_input<-data.frame(PCA$x[positions,1:2])
space4 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 301-400")

# make morpho space for 401-500
positions <- which(rownames(PCA$x) %in% depth.401TO500$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space5_input<-data.frame(PCA$x[positions,1:2])
space5 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 401-500")

# make morpho space for 501-600
positions <- which(rownames(PCA$x) %in% depth.501TO600$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space6_input<-data.frame(PCA$x[positions,1:2])
space6 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 501-600")

# make morpho space for 601-700
positions <- which(rownames(PCA$x) %in% depth.601TO700$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space7_input<-data.frame(PCA$x[positions,1:2])
space7 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 601-700")

# make morpho space for 701-800
positions <- which(rownames(PCA$x) %in% depth.701TO800$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space8_input<-data.frame(PCA$x[positions,1:2])
space8 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 701-800")

# create a data frame of the plots with the phytools components
space1_df <-phytools2ggplot(space1)
space2_df <-phytools2ggplot(space2)
space3_df <-phytools2ggplot(space3)
space4_df <-phytools2ggplot(space4)
space5_df <-phytools2ggplot(space5)
space6_df <-phytools2ggplot(space6)
space7_df <-phytools2ggplot(space7)
space8_df <-phytools2ggplot(space8)

# create a list of the plots
plots1 <- list(
  ggplot(data = space1_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space1_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 0-100")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space2_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space2_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 101-200")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space3_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space3_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 201-300")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space4_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space4_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 301-400")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space5_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space5_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 401-500")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space6_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space6_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 501-600")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space7_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space7_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 601-700")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    ,
  
  ggplot(data = space8_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space8_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 701-800")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
    
)
# plot the list of plots as a facet plot
dev.new(width=5, height=2, unit="in")
grid.arrange(grobs = plots1, ncol = 2, respect = TRUE)
#-------------------------------------------------------------------------------

# make facet plot for all morpho spaces by 200s
#-------------------------------------------------------------------------------
# make morpho space for 0-200 
positions <- which(rownames(PCA$x) %in% depth.0TO200$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space9_input<-data.frame(PCA$x[positions,1:2])
space9 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main="Depth: 0-200")

# make morpho space for 201-400
positions <- which(rownames(PCA$x) %in% depth.201TO400$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space10_input<-data.frame(PCA$x[positions,1:2])
space10 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 201-400")

# make morpho space for 401-600 
positions <- which(rownames(PCA$x) %in% depth.401TO600$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space11_input<-data.frame(PCA$x[positions,1:2])
space11 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label="radial",xlab="PCA1",ylab="PCA2")
title(main="Depth: 401-600")

# make morpho space for 601-800
positions <- which(rownames(PCA$x) %in% depth.601TO800$Taxon)
tipsTOdrop <- c(rownames(PCA$x))
diff <- setdiff(tipsTOdrop,rownames(PCA$x[positions,1:2]))
temp.tree <- drop.tip(new.tree,diff)
space12_input<-data.frame(PCA$x[positions,1:2])
space12 <- phylomorphospace(temp.tree,PCA$x[positions,1:2],label=FALSE,xlab="PCA1",ylab="PCA2")
title(main="Depth: 601-800")

# create a dataframe of the plots with the phytools components
space9_df <-phytools2ggplot(space9)
space10_df <-phytools2ggplot(space10)
space11_df <-phytools2ggplot(space11)
space12_df <-phytools2ggplot(space12)

# create a list of the plots
plots2 <- list(
  ggplot(data = space9_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space9_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 0-200")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  
  ggplot(data = space10_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space10_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 201-400")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  
  ggplot(data = space11_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space11_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 401-600")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1),
  
  ggplot(data = space12_input, aes(x = Comp1, y = Comp2)) +
    geom_point(size = 2, color = "darkblue") +
    geom_segment(data = space12_df, aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
                 linewidth = 0.5, color = "darkblue") +
    theme_classic() +
    ggtitle("Depth: 601-800")
    + xlim(-0.2, 0.1) + ylim(-0.1,0.1)
)

# plot the list of plots as a facet plot
dev.new(width=5, height=2, unit="in")
grid.arrange(grobs = plots2, ncol = 2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT4 = p4_depthRangeVsLengthPlot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# depth range vs length = rain cloud plot (pirate plot maybe)
#-------------------------------------------------------------------------------
# create new df based on the species from phylomorphospace
depth.data.filtered <- rbind(depth.0TO100,depth.101TO200,depth.201TO300,depth.301TO400,depth.401TO500,depth.501TO600,depth.601TO700,depth.701TO800)

# Create a vector to store the group identities
groups <- character(length(depth.data.filtered$Depth))

# Assign group identities based on data values
groups[depth.data.filtered$Depth >= 0 & depth.data.filtered$Depth <= 100] <- "0to100"
groups[depth.data.filtered$Depth >= 101 & depth.data.filtered$Depth <= 200] <- "101to200"
groups[depth.data.filtered$Depth >= 201 & depth.data.filtered$Depth <= 300] <- "201to300"
groups[depth.data.filtered$Depth >= 301 & depth.data.filtered$Depth <= 400] <- "301to400"
groups[depth.data.filtered$Depth >= 401 & depth.data.filtered$Depth <= 500] <- "401to500"
groups[depth.data.filtered$Depth >= 501 & depth.data.filtered$Depth <= 600] <- "501to600"
groups[depth.data.filtered$Depth >= 601 & depth.data.filtered$Depth <= 700] <- "601to700"
groups[depth.data.filtered$Depth >= 701 & depth.data.filtered$Depth <= 800] <- "701to800"

# Add the group column to the data frame
depth.data.filtered$group <- groups

#raincloud plot
dev.new(width=5, height=2, unit="in")
ggplot(depth.data.filtered, aes(y = Length, x = group, fill = group)) +
  stat_halfeye(
    adjust = 0.5,
    justification = -0.2,
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.color = NA,
    alpha = 0.5
  ) +
  geom_jitter(width = .05, alpha = .3) +
  scale_fill_tq() +
  theme_tq() +
  labs(
    title = "RainCloud Plot",
    x = "Depth",
    y = "Lengths",
    fill = "Depth"
  ) +
  coord_flip()
#-------------------------------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT5 = p5_iqrLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# IQR length vs depth range = regular plot (lollipop)
#-------------------------------------------------------------------------------
calc.inter.quartile <- function(df){
  df <- df %>% filter(!Length=="NA")
  return(IQR(df$Length))
}

# apply function to all the depth ranges by 100s
list.of.depth.data.100s <- list(depth.0TO100,depth.101TO200,depth.201TO300,depth.301TO400,depth.401TO500,depth.501TO600,depth.601TO700,depth.701TO800)
iqr.100s <- lapply(list.of.depth.data.100s,calc.inter.quartile)
print(iqr.100s)
iqr.100s <- data.frame(iqr.100s)
colnames(iqr.100s) <- list("0to100","101to200","201to300","301to400","401to500","501to600","601to700","701to800")
iqr.100s <- t(iqr.100s)
iqr.100s <- data.frame(iqr.100s)
iqr.100s <- tibble::rownames_to_column(iqr.100s, "x")
colnames(iqr.100s)[2] <- "y"

# plot
dev.new(width=5, height=2, unit="in")
ggplot(iqr.100s,aes(x = x , y = y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y)) +
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 

list.of.depth.data.200s <- list(depth.0TO200,depth.201TO400,depth.401TO600,depth.601TO800)
iqr.200s <- lapply(list.of.depth.data.200s,calc.inter.quartile)
print(iqr.200s)
iqr.200s <- data.frame(iqr.200s)
colnames(iqr.200s) <- list("0to200","201to400","401to600","601to800")
iqr.200s <- t(iqr.200s)
iqr.200s <- data.frame(iqr.200s)
iqr.200s <- tibble::rownames_to_column(iqr.200s, "x")
colnames(iqr.200s)[2] <- "y"

dev.new(width=5, height=2, unit="in")
ggplot(iqr.200s,aes(x = x , y = y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y)) +
  geom_point( size=5, color="blue", fill=alpha("green", 0.3), alpha=0.7, shape=21, stroke=2) 
#-------------------------------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT6 = p6_varianceLengthVsDepthRange
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Variance length vs depth range = regular plot (lollipop)
#-------------------------------------------------------------------------------
# Calculate coefficient of variance ( sd/mean ) same as IQR graph 
calc.variance <- function(df){
  df <- df %>% filter(!Length=="NA")
  return(sd(df$Length)/mean(df$Length))
}

# apply function to all the depth ranges by 100s
list.of.depth.data <- list(depth.0TO100,depth.101TO200,depth.201TO300,depth.301TO400,depth.401TO500,depth.501TO600,depth.601TO700,depth.701TO800)
var.100s <- lapply(list.of.depth.data,calc.variance)
print(var.100s)
var.100s <- data.frame(var.100s)
colnames(var.100s) <- list("0to100","101to200","201to300","301to400","401to500","501to600","601to700","701to800")
var.100s <- t(var.100s)
var.100s <- data.frame(var.100s)
var.100s <- tibble::rownames_to_column(var.100s, "x")
colnames(var.100s)[2] <- "y"

# plot
dev.new(width=5, height=2, unit="in")
ggplot(var.100s,aes(x = x , y = y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y)) +
  geom_point( size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 


list.of.depth.data.200s <- list(depth.0TO200,depth.201TO400,depth.401TO600,depth.601TO800)
var.200s <- lapply(list.of.depth.data.200s,calc.variance)
print(var.200s)
var.200s <- data.frame(var.200s)
colnames(var.200s) <- list("0to200","201to400","401to600","601to800")
var.200s <- t(var.200s)
var.200s <- data.frame(var.200s)
var.200s <- tibble::rownames_to_column(var.200s, "x")
colnames(var.200s)[2] <- "y"

dev.new(width=5, height=2, unit="in")
ggplot(var.200s,aes(x = x , y = y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y)) +
  geom_point( size=5, color="purple", fill=alpha("pink", 0.3), alpha=0.7, shape=21, stroke=2) 
#-------------------------------------------------------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT7 = p7_bodyLengthVsComps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# body length vs comp1 = phylomorphospace
#-------------------------------------------------------------------------------
change.name.fishBase <- function(data){
  for (i in 1:length(data$Taxon)){
    for (j in 1:length(fishBase.nameChange$Old)){
      if (data$Taxon[i] == fishBase.nameChange$Old[j]){
        data$Taxon[i]<-fishBase.nameChange$New[j]
      }
    }
  }
  # use FishBase to find valid names and make a new file with updated names
  species_names <- unique(data$Taxon)
  lengths <- c()
  for (i in species_names){
    dat <- species(i, fields=c("Length"))
    lengths <- c(lengths,dat)
  } 
  names(lengths) <- species_names
  if ("Nototheniops nudifrons" %in% names(lengths) ){
    lengths$`Nototheniops nudifrons`[1] <- 19.5
  }
  if ("Trematomus pennelii" %in% names(lengths) ){
    lengths$`Trematomus pennelii`[1] <- 29.7
  }
  if ("Trematomus loenbergii" %in% names(lengths) ){
    lengths$`Trematomus loenbergii`[1] <- 30.0
  }
  if ("Artedidraco skottsbergi" %in% names(lengths) ){
    lengths$`Artedidraco skottsbergi`[1] <- 12.0
  }
  if ("Cryodraco antarcticus" %in% names(lengths) ){
    lengths$`Cryodraco antarcticus`[1] <- 48.7
  }
  return(lengths)
}

change.name.PCA <- function(data){
  for (i in 1:length(names(data))){
    for (j in 1:length(fishBase.nameChange$Old)){
      if (names(data)[i] == fishBase.nameChange$Old[j]){
        names(data)[i]<-fishBase.nameChange$New[j]
      }
    }
  }
  return(data)
}
pca.data <- PCA$x[,1:2]
size.data <- change.name.fishBase(depth.data.filtered)
fishB.names <- fishBase.nameChange


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fishB.names$New<-gsub("[.]", " ", fishB.names$New)

matching_indices <- match(names(size.data), fishB.names$New)
SizeData <- size.data[matching_indices]

df <- do.call(rbind, SizeData)
SizeData.df <- as.data.frame(df)

fishB.names$SizeData<-SizeData.df$V1

names(fishB.names)<-c("New","X","Size")

pca.df <- data.frame(pca.data)
pca.df <- tibble::rownames_to_column(pca.df, "X") 
BSData<-fishB.names[,2:3]

merged_data<-merge(pca.df, BSData, by = "X", all.x = TRUE)
merged_data<-na.omit(merged_data)
rownames(merged_data) <- merged_data$X
merged_data[,4] <- log(merged_data[,4])
matrix_merged_data <- data.matrix(merged_data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rename the tree to match data
diff.tips.fishBase <- setdiff(new.tree$tip.label,merged_data$X)
fishBase.tree <- drop.tip(new.tree,diff.tips.fishBase)
print("NEW PHYLOGENETIC TREE COMPLETE")
dev.new(width=5, height=5, unit="in")
plot(fishBase.tree)

dev.new(width=5, height=5, unit="in")
phylomorphospace(fishBase.tree,merged_data[,c(2,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp1",font.main=3)

dev.new(width=5, height=5, unit="in")
phylomorphospace(fishBase.tree,merged_data[,c(3,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp2",font.main=3)
#-------------------------------------------------------------------------------


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT8 = p8_speciesVsDepth+speciesVsLength
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# species vs depth = regular plot (violin+dots) (last 4 slides in presentation)
#-------------------------------------------------------------------------------
dev.new()
ggplot(depth.data,aes(x= Depth,y=Taxon, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Depth", option = "C") +
  labs(title = "Depth of species") +
  geom_jitter(width = .05, alpha = .3)

#-------------------------------------------------------------------------------

# species vs length = regular plot (violin+dots) (last 4 slides in presentation)
#-------------------------------------------------------------------------------
depth.data <- na.omit(depth.data)
outliers <- c()
for (k in 1:nrow(depth.data)) {
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
  labs(title = "Length of species") +
  geom_jitter(width = .05, alpha = .3)

lepi_squa.df <- which(no.outlier.depth.data$Taxon=="lepi_squa-trematominae")
lepi.location <- which(fishB.names$X == "lepi_squa-trematominae")
size <- fishB.names[lepi.location,]$Size
lepi_squa.df <- no.outlier.depth.data[lepi_squa.df,]
new.df <- which(lepi_squa.df$Length>size)
lepi_squa.df <- lepi_squa.df[new.df,]

# check why this does not work
dev.new()
ggplot(depth.data.filtered,aes(x= Length,y=Taxon, fill = group)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Length", option = "C") +
  labs(title = "Length of species") +
  geom_jitter(width = .05, alpha = .3)
#-------------------------------------------------------------------------------
