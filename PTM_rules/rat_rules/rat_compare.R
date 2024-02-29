library(dplyr)

sampledf3 <- read.table("rat_sample.txt", sep = '\t', col.names = c("m.z", "intensity"))
sampledf3 <- sampledf3 %>%
  mutate(across(c('m.z'), round, 1))
sampledf3$m.z_plus1 <- sampledf3$m.z + 0.1
sampledf3$m.zminus1 <- sampledf3$m.z - 0.1

rat_table_df <- read.csv("results.csv", sep = ",")
rat_table_df <- rat_table_df %>%
  mutate(across(c(mass1), round, 1))

mass_list = c()
count = 0
for (mass in rat_table_df$mass1) {
  for (row in 1:nrow(sampledf3)) {
    if (mass <= sampledf3[row, "m.z_plus1"] & mass >= sampledf3[row,"m.zminus1"]) {
      mass_list <- append(mass_list, sampledf3[row, "m.z"])
      #print(mass)
      #print(sampledf[row, "m.z_plus1"])
      #print(sampledf[row,"m.zminus1"])
      count = count + 1
    } 
  }
}
print(count)
print(length(unique(mass_list)))