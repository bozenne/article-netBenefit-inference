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

data.figure1 <- readRDS(file = file.path(path.results,"dataTiming-SIMULATION_H0-1TTE.rds"))
## data.figure1[,.N,by = c("n","method", "type")]
## based on only one file i.e. a subset of the simulations (100 simulations)

## * Create figure
figure1 <- ggTiming(data.figure1, type.data = "processed", plot = FALSE) 

## * Export
ggsave(figure1$plot, filename = file.path("figures-article","figure1.pdf"), width = 10, height = 9, device = cairo_pdf)

## figure1$data[n==1000,mean(time),by = c("method","type")]
##    method               type        V1
## 1:  Gehan              ratio  1.095203
## 2:  Peron              ratio 18.363595
## 3:  Gehan estimate with s.e.  0.126220
## 4:  Peron estimate with s.e.  7.846810
## 5:  Gehan           estimate  0.115260
## 6:  Peron           estimate  0.427330
