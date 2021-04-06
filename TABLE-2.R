## * path
## path <- "P:/Cluster/GPC/Article-inference-Ustatistic-Rao"
## setwd(path)
path.results <- "Results"
path.tables <- "tables-article"


## * libraries
library(data.table)
library(xtable) ## neededed for display
source("FCT-gg.R")

## * Load results

data.table2 <- readRDS(file = file.path(path.results,"dataTable-SIMULATION_H1-1TTE.rds"))

## * Create table
table2 <- createTable(data.table2, type.data = "processed", print = FALSE,
                      digits = 3, by = "threshold") 

## * export
fileConn <- file(file.path(path.tables,"table2.txt"))
writeLines(table2$table, fileConn)
close(fileConn)
