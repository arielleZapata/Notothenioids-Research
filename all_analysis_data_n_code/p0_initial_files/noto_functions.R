# make morpho spaces + calculate convex hulls function
calcConvex.phylomorpho.chull <- function(df,title,new.tree){
  rownames(PCA) <- PCA$X
  positions <- which(PCA$X %in% df$Taxon)
  tipsTOdrop <- PCA$X
  diff <- setdiff(tipsTOdrop,rownames(PCA[positions,2:3]))
  td <- treedata(new.tree,PCA[positions,2:3])  
  temp.tree <- drop.tip(new.tree,diff)
  dev.new()
  space<-phylomorphospace(td$phy,td$data,label=FALSE,xlab="PCA1",ylab="PCA2")
  title(main=title)
  pc <- as_tibble(PCA[positions,2:3])
  chull<-chull(pc)
  convex_hull <- pc %>% slice(chull)
  area <- as.matrix(convex_hull) %>% areapl
  output<-list(space,area)
  return(output)
}
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
# filter out NAs and calculate IQR
calc.inter.quartile <- function(df){
  df <- df %>% dplyr::filter(!Length=="NA")
  return(IQR(df$Length))
}
# filter out NAs and calculate variance
calc.variance <- function(df){
  df <- df %>% dplyr::filter(!Length=="NA")
  return(sd(df$Length)/mean(df$Length))
}
# change names in data to match fishBase names
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
# change names in PCA to match fishBase names
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
# make morpho space for depth range
make.morpho.depth <- function(PCA,Range=depth.0TO100,tree=new.tree,title="Depth: 0-100",...){
positions <- which(PCA$X %in% depth.0TO100$Taxon)
PCAsubset <- PCA[positions,]
tipsTOdrop <- c(PCA$X)
diff <- setdiff(tipsTOdrop,PCA$X[positions])
temp.tree <- drop.tip(new.tree,diff)
diff2 <- temp.tree$tip.label[which(!temp.tree$tip.label%in%PCA$X)]
temp.tree <- drop.tip(temp.tree,diff2)
space1_input <- data.frame(cbind(PCAsubset$Comp1,PCAsubset$Comp2))
rownames(space1_input) <- PCAsubset$X
space1 <- phylomorphospace(temp.tree,space1_input,label="radial",xlab="PCA1",ylab="PCA2")
title(title)
return(list(space1_input,space1))
}