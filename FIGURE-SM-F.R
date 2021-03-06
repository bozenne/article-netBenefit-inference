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
data.figureSMF.timing <- readRDS(file = file.path(path.results,"dataTiming-SIMULATION_H0-mE.rds"))

data.figureSMF.H0a <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-1TTE.rds"))
data.figureSMF.H0b <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-mE.rds"))

data.figureSMF.H1a <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H1-1TTE.rds"))
data.figureSMF.H1b <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H1-mE.rds"))

data.figureSMF.nmH0a <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-1TTE-nm.rds"))
data.figureSMF.nmH0b <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H0-mE-nm.rds"))

data.figureSMF.nmH1a <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H1-1TTE-nm.rds"))
data.figureSMF.nmH1b <- readRDS(file = file.path(path.results,"dataCoverage-SIMULATION_H1-mE-nm.rds"))

## * Create figure
figureSMF.timing <- ggTiming(data.figureSMF.timing, type.data = "processed", plot = FALSE) 

figureSMF.H0a <- ggCoverage(data.figureSMF.H0a, type.data = "processed", plot = FALSE) 
figureSMF.H0b <- ggCoverage(data.figureSMF.H0b, type.data = "processed", plot = FALSE) 
figureSMF.H0 <- ggarrange(figureSMF.H0a$plot + coord_cartesian(ylim = c(0.9,1)),
                          figureSMF.H0b$plot + coord_cartesian(ylim = c(0.9,1)),
                          ncol = 1,
                          common.legend = TRUE, legend = "bottom")

figureSMF.H1a <- ggCoverage(data.figureSMF.H1a, type.data = "processed", plot = FALSE) 
figureSMF.H1b <- ggCoverage(data.figureSMF.H1b, type.data = "processed", plot = FALSE) 
figureSMF.H1 <- ggarrange(figureSMF.H1a$plot + coord_cartesian(ylim = c(0.9,1)),
                          figureSMF.H1b$plot + coord_cartesian(ylim = c(0.9,1)),
                          ncol = 1,
                          common.legend = TRUE, legend = "bottom")

figureSMF.nmH0a <- ggCoverage(data.figureSMF.nmH0a, type.data = "processed", plot = FALSE) 
figureSMF.nmH0b <- ggCoverage(data.figureSMF.nmH0b, type.data = "processed", plot = FALSE) 
figureSMF.nmH0 <- ggarrange(figureSMF.nmH0a$plot + coord_cartesian(ylim = c(0.9,1)),
                            figureSMF.nmH0b$plot + coord_cartesian(ylim = c(0.9,1)),
                            ncol = 1,
                            common.legend = TRUE, legend = "bottom")

figureSMF.nmH1a <- ggCoverage(data.figureSMF.nmH1a, type.data = "processed", plot = FALSE) 
figureSMF.nmH1b <- ggCoverage(data.figureSMF.nmH1b, type.data = "processed", plot = FALSE) 
figureSMF.nmH1 <- ggarrange(figureSMF.nmH1a$plot + coord_cartesian(ylim = c(0.9,1)),
                            figureSMF.nmH1b$plot + coord_cartesian(ylim = c(0.9,1)),
                            ncol = 1,
                            common.legend = TRUE, legend = "bottom")


## * Export
ggsave(figureSMF.timing$plot, filename = file.path("figures-article","figureSMF-timing.pdf"), width = 10, height = 9, device = cairo_pdf)
## figureSMF.timing$data[n==1000, mean(time),by = c("scoring.rule","type")]
##           scoring.rule               type        V1
## 1: Gehan's scoring rule              ratio  1.110963
## 2: Peron's scoring rule              ratio 20.863980
## 3: Gehan's scoring rule estimate with s.e.  0.269020
## 4: Peron's scoring rule estimate with s.e. 14.366780
## 5: Gehan's scoring rule           estimate  0.242200
## 6: Peron's scoring rule           estimate  0.688660

ggsave(figureSMF.H0, filename = file.path("figures-article","figureSMF-H0.pdf"), width = 10, height = 9, device = cairo_pdf)
ggsave(figureSMF.H1, filename = file.path("figures-article","figureSMF-H1.pdf"), width = 10, height = 9, device = cairo_pdf)
ggsave(figureSMF.nmH0, filename = file.path("figures-article","figureSMF-nmH0.pdf"), width = 10, height = 9, device = cairo_pdf)
ggsave(figureSMF.nmH1, filename = file.path("figures-article","figureSMF-nmH1.pdf"), width = 10, height = 9, device = cairo_pdf)


