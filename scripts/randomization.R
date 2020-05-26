raw_data <- readRDS("data/raw_data.rds")

set.seed(20191219)

random_data <- data.frame(conc = raw_data$conc)

for (i in unique(raw_data$conc)) {
  for (j in colnames(raw_data)[-1]) {
    random_data[random_data$conc == i, j] <- sample(raw_data[raw_data$conc == i, j])
  }
}

write.csv(random_data, "data/random_data.csv", row.names = FALSE)
saveRDS(random_data, "data/random_data.rds")
