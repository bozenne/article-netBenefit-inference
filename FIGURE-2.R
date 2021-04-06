## * path
## path <- "P:/Cluster/GPC/Article-inference-Ustatistic-Rao"
## setwd(path)
path.results <- "Results"
path.figures <- "figures-article"


## * libraries
library(data.table)
library(ggplot2)
library(ggpubr) ## neededed for graphical display
library(ggthemes) ## neededed for graphical display
library(xtable) ## neededed for display
source("FCT-gg.R")

## * Load results

data.figure2 <- readRDS(file = file.path(path.results,"dataTiming-SIMULATION_H0-1TTE.rds"))
## data.figure2[,.N,by = c("n","method", "type")]
## based on only one file i.e. a subset of the simulations (100 simulations)

## * Create figure
figure2 <- ggTiming(data.figure2, type.data = "processed", plot = FALSE) 

## * Export
ggsave(figure2$plot, filename = file.path("figures-article","figure2.pdf"), width = 10, height = 9, device = cairo_pdf)
