#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DDAP6 = DDAP6_bodyLengthVsComps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(rfishbase)
library(phytools)

# inputs
source("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
PCA <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/PCAcomps_OUTPUTS.csv") # csv of PCA data
depth.data.filtered <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p9_csvAllLocations/SouthOrkneyIs/DDSOI1_pholomorpho+depthPlots/DDSOI.depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data
fishBase.nameChange <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/fishBase_names.csv") # csv of name replacement for retrieving the fishbase data
new.tree <- read.tree("~/Notothenioids-Research/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/new.tree.tre") # import trimmed phylogenetic tree
csv.names <- read.csv("~/Notothenioids-Research/all_analysis_data_n_code/p0_initial_files/updated_names.csv",header=FALSE) # csv of name replacement for tree

# match PCA names to fishBase names
colnames(depth.data.filtered)[colnames(depth.data.filtered) == "Taxon"] <- "Species"
pca.data <- PCA
size.data <- change.name.fishBase(depth.data.filtered)
fishB.names <- fishBase.nameChange

fishB.names$New<-gsub("[.]", " ", fishB.names$New)
matching_indices <- match(names(size.data), fishB.names$New)
SizeData <- size.data[matching_indices]
df <- do.call(rbind, SizeData)
SizeData.df <- as.data.frame(df)
names(SizeData.df) <- "SizeData"
SizeData.df$New <- rownames(SizeData.df)

mergedSizeData <- merge(SizeData.df,fishB.names,by="New",all.x=TRUE)
pca.df <- data.frame(pca.data)
names(pca.df) <- c("Old","Comp1","Comp2")

# merge the data at column X
merged_data <- merge(pca.df, mergedSizeData, by = "Old", all.x = TRUE)
merged_data <- na.omit(merged_data)
rownames(merged_data) <- merged_data$Old
merged_data$LogSize <- log(merged_data$SizeData)
matrix_merged_data <- data.matrix(merged_data)

# rename the tree to match data
tips <- c()
for(i in 1:nrow(csv.names)) {
  if (any(csv.names[i, 2]=="")){
    tips <- c(tips,csv.names[i,1])
  }
  else{
    new.tree$tip.label[new.tree$tip.label==csv.names[i, 1]] <- csv.names[i, 2]
  }
}
diff.tips.fishBase <- setdiff(new.tree$tip.label,merged_data$Old)
fishBase.tree <- drop.tip(new.tree,diff.tips.fishBase)
dev.new()
plot(fishBase.tree)

colnames(matrix_merged_data) <- names(merged_data)

# display the new phylomorphospaces
dev.new()
output1 <- phylomorphospace(fishBase.tree,matrix_merged_data[,c("Comp1","LogSize")],label=FALSE)
title(main="Phylomorphospace Plot of Comp1",font.main=3)

output2 <- dev.new(width=5, height=5, unit="in")
phylomorphospace(fishBase.tree,matrix_merged_data[,c("Comp2","LogSize")],label=FALSE)
title(main="Phylomorphospace Plot of Comp2",font.main=3)

## ----------outputs----------
### CSV
write.csv(merged_data,file="merged_data.csv")

### PDF
pdf(file= "DDSOI6.bodyLengthVsComps_OUTPUTS.pdf")

phylomorphospace(fishBase.tree,matrix_merged_data[,c("Comp1","LogSize")],label=FALSE)
title(main="Phylomorphospace Plot of Comp1",font.main=3)

phylomorphospace(fishBase.tree,matrix_merged_data[,c("Comp2","LogSize")],label=FALSE)
title(main="Phylomorphospace Plot of Comp2",font.main=3)

plot(fishBase.tree)

dev.off() 
