#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PT7 = p7_bodyLengthVsComps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# libraries
library(rfishbase)
library(phytools)

# inputs
source("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/noto_functions.R")
PCA <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p1_gpa_phylomorpho_all_data/PCAcomps_OUTPUTS.csv") # csv of PCA data
depth.data.filtered <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p2_pholomorpho+depthPlots/depth.data.filtered_OUTPUTS.csv") # csv of depth.data.filtered data
fishBase.nameChange <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/fishBase_names.csv") # csv of name replacement for retrieving the fishbase data
new.tree <- read.tree("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/notothenioid_timetree.tre") # import trimmed phylogenetic tree
csv.names <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/updated_names.csv",header=FALSE) # csv of name replacement for tree

# ADD FILTERED DATA FOR EACH LOCATION

# match PCA names to fishBase names
pca.data <- PCA
size.data <- change.name.fishBase(depth.data.filtered)
fishB.names <- fishBase.nameChange

fishB.names$New<-gsub("[.]", " ", fishB.names$New)
matching_indices <- match(names(size.data), fishB.names$New)
SizeData <- size.data[matching_indices]
df <- do.call(rbind, SizeData)
SizeData.df <- as.data.frame(df)
fishB.names$SizeData<-SizeData.df$V1
names(fishB.names)<-c("New","X","Size")
pca.df <- data.frame(pca.data)
BSData <- fishB.names[,2:3]

# merge the data at column X
merged_data <- merge(pca.df, BSData, by = "X", all.x = TRUE)
merged_data <- na.omit(merged_data)
rownames(merged_data) <- merged_data$X
merged_data[,4] <- log(merged_data[,4])
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
diff.tips.fishBase <- setdiff(new.tree$tip.label,merged_data$X)
fishBase.tree <- drop.tip(new.tree,diff.tips.fishBase)
dev.new()
plot(fishBase.tree)

# display the new phylomorphospaces
dev.new()
phylomorphospace(fishBase.tree,merged_data[,c(2,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp1",font.main=3)

dev.new(width=5, height=5, unit="in")
phylomorphospace(fishBase.tree,merged_data[,c(3,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp2",font.main=3)

## ----------outputs----------
### CSV
merged_data.revLog <- merged_data
merged_data.revLog$Size <- exp(merged_data.revLog$Size)
write.csv(merged_data.revLog,file="merged_data.revLog.csv")

### PDF
pdf(file= "bodyLengthVsComps_OUTPUTS.pdf")
phylomorphospace(fishBase.tree,merged_data[,c(2,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp1",font.main=3)

phylomorphospace(fishBase.tree,merged_data[,c(3,4)],label=FALSE)
title(main="Phylomorphospace Plot of Comp2",font.main=3)

plot(fishBase.tree)

dev.off() 