#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT1 = p1_gpa_and_phylomorpho_of_all_data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(geomorph)
library(phytools)

#set wd to folder with needed content
setwd("~/Notothenioids-Research/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data")

# make variable for file locations
tps.file <- "~/Notothenioids-Research/all_analysis_data_n_code/make_tps_files+combined_p_file/all_tps_file.tps" 

# read in data
fish <- read.tree("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/notothenioid_timetree.tre") # given phylo tree
csv.names <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/updated_names.csv",header=FALSE) # csv of name replacement for tree

# read the TPS file
data <- readland.tps(tps.file,specID = c("None", "ID", "imageID"), negNA = FALSE,readcurves = FALSE, warnmsg = FALSE)

# perform and display GPA on TPS file
data.gpa <- gpagen(data,print.progress = FALSE)
summary(data.gpa)
dev.new()
plotAllSpecimens(data.gpa$coords)
title("GPA Consensus")

# read in and display time tree
fish.phylo <- as.phylo(fish)
dev.new()
plot(fish.phylo)
title("Given Phylogenetic Tree for Notothenioids")

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
dev.new()
plot(new.tree)
title("Pruned Given Phylogenetic Tree for Notothenioids")

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
dev.new()
plot(PCA$x[,1:2])
title("PCA - Comp1 vs Comp2")

# plot the data in phylomorphospace
dev.new()
phylomorphospace(new.tree,PCA$x[,1:2],label=FALSE)
title(main="Phylomorphospace Plot of All Data",font.main=3)

## ----------outputs----------
# CSV
write.csv(PCA$x[,1:2], "PCAcomps_OUTPUTS.csv")

# TRE
write.tree(new.tree, file = "new.tree.tre")

# PDF
pdf(file= "gpaCons+phylotree+pcaComp_OUTPUTS.pdf" )
plotAllSpecimens(data.gpa$coords)
title("GPA Consensus")
plot(fish.phylo)
title("Given Phylogenetic Tree for Notothenioids")
plot(new.tree)
title("Pruned Given Phylogenetic Tree for Notothenioids")
plot(PCA$x[,1:2])
title("PCA - Comp1 vs Comp2")
phylomorphospace(new.tree,PCA$x[,1:2],label=FALSE)
title(main="Phylomorphospace Plot of All Data",font.main=3)
dev.off() 