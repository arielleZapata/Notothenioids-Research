# create a new csv for each location in Combined Antartic data
depth.data <- read.csv("~/Notothenioids_research/repository/all_analysis_data_n_code/p0_initial_files/CombinedAntarcticData.csv") # depth data for fish

# see how many and which locations are recorded
unique(depth.data$Region)

# subset the data sets by location
depth.data.South.Shetland.Is <- subset(depth.data, Region == 'South Shetland Is.')
unique(depth.data.South.Shetland.Is$Region)

depth.data.South.Orkney.Is.1 <- subset(depth.data, Region == 'South Orkney Is.')
depth.data.South.Orkney.Is.2 <- subset(depth.data, Region == 'South Orkney Island')
depth.data.South.Orkney.Is <- rbind(depth.data.South.Orkney.Is.1,depth.data.South.Orkney.Is.2)
unique(depth.data.South.Orkney.Is$Region)

depth.data.Antarctic.Peninsula <- subset(depth.data, Region == 'Antarctic Peninsula')
unique(depth.data.Antarctic.Peninsula$Region)

depth.data.Elephant.Island <- subset(depth.data, Region == 'Elephant Island')
unique(depth.data.Elephant.Island$Region)

# save all the df to csv files
write.csv(depth.data.South.Shetland.Is,"~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.South.Shetland.Is.csv")
write.csv(depth.data.South.Orkney.Is,"~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.South.Orkney.Is.csv")
write.csv(depth.data.Antarctic.Peninsula,"~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.Antarctic.Peninsula.csv")
write.csv(depth.data.Elephant.Island,"~/Notothenioids_research/repository/all_analysis_data_n_code/p9_csvAllLocations/depth.data.Elephant.Island.csv")