##### Prepare simulation results for plot

## Function to calculate mean and CI of estimation errors
get_mean_CI <- function(setting, type, nrepl) {
  if (setting == "1") {
    filename <- paste0("S1_n150ones_estErr", type, ".txt")
    xVar <- -3:5
  }
  if (setting == "2") {
    filename <- paste0("S2_varyn_estErr", type, ".txt")
    xVar <- seq(from = 100, by = 50, to = 300)
  }
  if (setting == "3") {
    filename <- paste0("S3_realAAYeo_estErr", type, ".txt")
    xVar <- -3:5
  }
  err <- read.table(file = filename, sep = "\t")
  err_mean <- colMeans(err)
  err_sd <- apply(err, 2, sd)
  err_lower <- err_mean - err_sd * 1.96 / sqrt(nrepl)
  err_upper <- err_mean + err_sd * 1.96 / sqrt(nrepl)
  return(data.frame(xVar = xVar, type = type, 
                    err_mean = err_mean, err_lower = err_lower, err_upper = err_upper))
}

### Setting 1
nrepl <- 100

results_setting1 <- NULL
for (t in c("S", "E", "N", "L", "R", "C")) {
  new_results <- get_mean_CI("1", t, nrepl)
  results_setting1 <- rbind(results_setting1, new_results)
}
colnames(results_setting1)[1] <- "k"


### Setting 2
results_setting2 <- NULL
for (t in c("S", "E", "N", "L", "R", "C")) {
  new_results <- get_mean_CI("2", t, nrepl)
  results_setting2 <- rbind(results_setting2, new_results)
}
colnames(results_setting2)[1] <- "n"


### Setting 3
results_setting3 <- NULL
for (t in c("S", "E", "N", "L", "R", "C")) {
  new_results <- get_mean_CI("3", t, nrepl)
  results_setting3 <- rbind(results_setting3, new_results)
}
colnames(results_setting3)[1] <- "k"


save(results_setting1, results_setting2, results_setting3,
     file = 'SimulationResults.Rdata')
