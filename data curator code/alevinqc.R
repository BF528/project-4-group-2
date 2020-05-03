featuredump = read.table("combinedSRR_output/alevin/featureDump.txt", sep = "\t", header = TRUE)
mean_maprate = mean(featuredump$MappingRate)
std_maprate = sd(featuredump$MappingRate)

mean_maprate
std_maprate

hist(featuredump$MappingRate, main = "mapping rates")