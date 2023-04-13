# load in all libraries needed
library(geomorph)
library(phytools)
library(dplyr)
library(magrittr)
library(ggplot2)
library(splancs)

# make variable for file locations
tps.file <- "~/Notothenioids_research/all_analysis_data_n_code/all_tps_file.tps" 

# read in data
fish <- read.tree("~/Notothenioids_research/all_analysis_data_n_code/notothenioid_timetree.tre") # given phylo tree
csv.names <- read.csv("~/Notothenioids_research/all_analysis_data_n_code/updated_names.csv",header=FALSE) # csv of name replacement for tree
depth.data <- read.csv("~/Notothenioids_research/all_analysis_data_n_code/CombinedAntarcticData.csv") # depth data for fish
depthNameChange <- read.csv("~/Notothenioids_research/all_analysis_data_n_code/FIXEDdepthNames.csv") # csv of name replacement for depth data


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

# names of the plot labels by the 100s
names <- c("Depth: 0-100","Depth: 101-200","Depth: 201-300","Depth: 301-400","Depth: 401-500","Depth: 501-600",
         "Depth: 601-700","Depth: 701-800")

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
plot(area.val,main="Depth Range vs CHull Area by the 100s",xlab="Depth Range (ft)",ylab="CHull Area")

# run function on 0-800 by 200s
b1 <- calcConvex.phylomorpho.chull(depth.0TO200,"Depth: 0-200")
b2 <- calcConvex.phylomorpho.chull(depth.201TO400,"Depth: 201-400")
b3 <- calcConvex.phylomorpho.chull(depth.401TO600,"Depth: 401-600")
b4 <- calcConvex.phylomorpho.chull(depth.601TO800,"Depth: 601-800")

# calculate and plot the areas for the 200s
area.val <- c(b1[[2]],b2[[2]],b3[[2]],b4[[2]])
dev.new(width=5, height=2, unit="in")
plot(area.val,main="Depth Range vs CHull Area by the 200s",xlab="Depth Range (ft)",ylab="CHull Area")

# Body size to depth
# like area calc interquartile range of the body size 
# facet plot all morpho spaces