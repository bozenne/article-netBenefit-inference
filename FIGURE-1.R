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
source("FCT-gg.R")

## * Load results

data.figure1a <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-1TTE.rds"))
data.figure1b <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-mE.rds"))
## data.figure1[,.N,by = c("n","method", "type")]
## based on only one file i.e. a subset of the simulations (100 simulations)

## * Create figure
figure1a <- ggCoverage(data.figure1a, type.data = "processed", plot = FALSE) 
figure1b <- ggCoverage(data.figure1b, type.data = "processed", plot = FALSE) 

figure1 <- ggarrange(figure1a$plot + coord_cartesian(ylim = c(0.9,1)),
                     figure1b$plot + coord_cartesian(ylim = c(0.9,1)),
                     ncol = 1,
                     common.legend = TRUE, legend = "bottom")

## * Export
ggsave(figure1, filename = file.path("figures-article","figure1.pdf"), width = 10, height = 9, device = cairo_pdf)
