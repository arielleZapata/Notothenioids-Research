# create folder to hold all tps files
setwd("~/Notothenioids_research/all_analysis_data_n_code")
dir.create("tps_files")

# find consensus files folder "official_consensus_files" and set working directory
main_dir <- "~/Notothenioids_research/all_analysis_data_n_code/official_consensus_files"
tps <- "~/Notothenioids_research/all_analysis_data_n_code/tps_files"
setwd(main_dir)

# get all folder names from the consensus files
folder.names <- list.files(main_dir)
folder.names

# loop through all the folders and find all the files
# the loop will combine them all to one tps file used later for GPA
setwd(tps)
sink("~/Notothenioids_research/all_analysis_data_n_code/tps_files/all_tps_file.tps")
# get into files needed
setwd(main_dir)
for (i in 1:length(folder.names)){
  path <- paste(main_dir,folder.names[i],sep="/")
  setwd(path)
  folder.names2 <- list.files(path)
  tps.temp <- paste(tps,folder.names[i],sep="/")
  tps.con <- paste(tps.temp,".tps",sep="")
  for (j in 1:length(folder.names2)){
    path2 <- paste(path,folder.names2[j],sep="/")
    lines <- readLines(path2)
    for (line in lines){
      cat(line)
      cat("\n")
    }
  }
}
sink()