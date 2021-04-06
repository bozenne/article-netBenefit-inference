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

data.table1 <- readRDS(file = file.path(path.results,"dataTable-SIMULATION_H0-1TTE.rds"))

## * Create table
table1 <- createTable(data.table1, type.data = "processed", print = FALSE,
                      digits = 3, by = "threshold") 

## * export
fileConn <- file(file.path(path.tables,"table1.txt"))
writeLines(table1$table, fileConn)
close(fileConn)
