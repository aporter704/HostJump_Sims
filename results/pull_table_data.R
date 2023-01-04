setwd("~/Dropbox/projects_WORKING/speciesJumpSimulations/HostJump_Sims/results/")
library(RSQLite)

conn <- dbConnect(SQLite(), 'MPXV_results.db')

dbGetInfo(conn)

tables_list <- dbListTables(conn)

# Here I will compare the number of jumps between all trees
tables_common <- grep('_common_', tables_list, value = T)
tables_common

jumps_results <- matrix(NA, length(tables_common), 9)

for(i in 1:length(tables_common)){
  tables_temp <- tables_common[i]
  command_query <- paste0("SELECT jumps, tree_name FROM ", tables_temp)
  query <- dbGetQuery(conn, command_query)
  if(nrow(query) == 9){
    jumps_results[i, ] <- query[, 1]
  }
}

colnames(jumps_results) <- query[, 2]
jumps_results <- as.data.frame(jumps_results)

dbDisconnect(conn)

colnames(jumps_results)
dim(jumps_results)

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(jumps_results$simulated_tree, jumps_results$simulated_non_painted_tree, 
     pch = 20, col = rgb(0, 0, 0, 0.1), 
     ylim = range(jumps_results$simulated_high_sampling_tree, na.rm = T),
     xlim = range(jumps_results$simulated_high_sampling_tree, na.rm = T))
abline(0, 1)
plot(jumps_results$simulated_tree, jumps_results$simulated_high_sampling_tree, 
     pch = 20, col = rgb(0, 0, 0, 0.1), 
     ylim = range(jumps_results$simulated_high_sampling_tree, na.rm = T),
     xlim = range(jumps_results$simulated_high_sampling_tree, na.rm = T))
abline(0, 1)
plot(jumps_results$simulated_tree, jumps_results$simulated_low_sampling_tree, 
     pch = 20, col = rgb(0, 0, 0, 0.1),
    ylim = range(jumps_results$simulated_high_sampling_tree, na.rm = T),
    xlim = range(jumps_results$simulated_high_sampling_tree, na.rm = T))
abline(0, 1)
plot(jumps_results$simulated_tree, jumps_results$simulated_opp_sampling_tree, 
    pch = 20, col = rgb(0, 0, 0, 0.1), 
    ylim = range(jumps_results$simulated_high_sampling_tree, na.rm = T),
    xlim = range(jumps_results$simulated_high_sampling_tree, na.rm = T))
abline(0, 1)

# Need to quantify error too. Perhaps as 