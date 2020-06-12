### BUILD_results.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 17 2018 (09:40) 
## Version: 
## Last-Updated: jun 12 2020 (16:08) 
##           By: Brice Ozenne
##     Update #: 124
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * path
library(data.table)
library(ggplot2)
path <- "P:/Cluster/GPC/Article-inference-Ustatistic"
setwd(path)
source("FCT-gg.R")

path.results <- file.path(path,"Results")

path.coverageH0_1TTE <- file.path(path.results,"CoverageH0-1TTE")
path.coverageH1_1TTE  <- file.path(path.results,"CoverageH1-1TTE")
path.coverageH0_mE <- file.path(path.results,"CoverageH0-mE")
path.coverageH1_mE  <- file.path(path.results,"CoverageH1-mE")


## * Coverage under the null (1 TTE)
cat("\n Coverage under the null (1 TTE) \n")
## ** load
dt.H0_1TTE <- butils::sinkDirectory(path.coverageH0_1TTE, string.keep = "tempo")
dt.H0_1TTE[, endpoint := NULL]

## ** display
ggTiming.H0_1TTE <- ggTiming(dt.H0_1TTE[threshold==0], file = 1)
ggBias.H0_1TTE  <- ggBias(dt.H0_1TTE[threshold==0.5])
ggSe.H0_1TTE  <- ggSe(dt.H0_1TTE[threshold==0])
ggCoverage.H0_1TTE  <- ggCoverage(dt.H0_1TTE[Hprojection==1 & threshold <= 0.5])
## ggCoverage.H0_1TTE$data
tableH0_1TTE  <- paste0(capture.output(
    createTable(dt.H0_1TTE[n %in% c(30,120,600) & Hprojection==1 & threshold <= 0.5],
                digits = 3)
), collapse = "\n")
cat(tableH0_1TTE)

## Hprojection effect
dt.H0_1TTE[threshold == 0.5 & method == "Gehan", .(sigma_empirical = sd(estimate),sigma_Ustat = mean(se)), by = c("Hprojection","n")]

## ** export
ggsave(ggTiming.H0_1TTE$plot, filename = file.path(path.results,"fig-timing-H0-1TTE.pdf"),
       width = 10, height = 9, device = cairo_pdf)
ggsave(ggBias.H0_1TTE$plot, filename = file.path(path.results,"fig-bias-H0-1TTE.pdf"),
       device = cairo_pdf)
ggsave(ggSe.H0_1TTE$plot, filename = file.path(path.results,"fig-se-H0-1TTE.pdf"),
       device = cairo_pdf)
ggsave(ggCoverage.H0_1TTE$plot + coord_cartesian(ylim = c(0.9,1)),
       filename = file.path(path.results,"fig-coverage-H0-1TTE.pdf"),
       height = 7, width = 10, device = cairo_pdf)
sink(file.path(path.results,"table-H0-1TTE.txt"))
cat(tableH0_1TTE)
sink()
saveRDS(dt.H0_1TTE, file = file.path(path.results,"dataCoverage.H0_1TTE.rds"))

## * Coverage under the alternative (1 TTE)
cat("\n Coverage under the alternative (1 TTE) \n")
## ** load
dt.H1_1TTE <- butils::sinkDirectory(path.coverageH1_1TTE, string.keep = "tempo")
dt.H1_1TTE[, endpoint := gsub("timeU","time",endpoint)]
## dt.H1_1TTE[n==400,.("all"= median(timeAll),"estimate"=mean(timeEstimate)),
## by = c("threshold","n","Hprojection","method")]
## dt.H1_1TTE <- readRDS(file = file.path(path.results,"dataCoverage.H1_1TTE.rds"))
dt.H1_1TTE[, endpoint := NULL]

## ** display
ggTiming.H1_1TTE <- ggTiming(dt.H1_1TTE[threshold==0], file = 1)
ggBias.H1_1TTE  <- ggBias(dt.H1_1TTE[threshold==0.5])
ggSe.H1_1TTE  <- ggSe(dt.H1_1TTE[threshold==0])
ggCoverage.H1_1TTE  <- ggCoverage(dt.H1_1TTE[Hprojection==1 & threshold <= 0.5])
tableH1_1TTE  <- paste0(capture.output(
    createTable(dt.H1_1TTE[n %in% c(30,120,600) & Hprojection==1 & threshold <= 0.5],
                digits = 3)
), collapse = "\n")
cat(tableH1_1TTE)
    
## ** export
ggsave(ggTiming.H1_1TTE$plot, filename = file.path(path.results,"fig-timing-H1-1TTE.pdf"),
       width = 10, height = 9, device = cairo_pdf)
ggsave(ggBias.H1_1TTE$plot, filename = file.path(path.results,"fig-bias-H1-1TTE.pdf"),
       device = cairo_pdf)
ggsave(ggSe.H1_1TTE$plot, filename = file.path(path.results,"fig-se-H1-1TTE.pdf"),
       device = cairo_pdf)
ggsave(ggCoverage.H1_1TTE$plot + coord_cartesian(ylim = c(0.9,1)),
       filename = file.path(path.results,"fig-coverage-H1-1TTE.pdf"),
       height = 7, width = 10, device = cairo_pdf)
sink(file.path(path.results,"table-H1-1TTE.txt"))
cat(tableH1_1TTE)
sink()
saveRDS(dt.H1_1TTE, file = file.path(path.results,"dataCoverage.H1_1TTE.rds"))


## * Coverage under the null (multiple TTE)
cat("\n Coverage under the null (multiple TTE) \n")

## ** load
dt.H0_mE <- butils::sinkDirectory(path.coverageH0_mE, string.keep = "tempo")
dt.H0_mE[, endpoint := gsub("timeU","time",endpoint)]
dt.H0_mE[, endpoint := factor(endpoint,
                              levels = c("time_0.5", "toxicity_0.5", "time_1e-12"),
                              labels = c("1 endpoint", "2 endpoints", "3 endpoints"))]

## ** display
ggTiming.H0_mE <- ggTiming(dt.H0_mE[endpoint=="3 endpoints"], file = 3)
ggBias.H0_mE  <- ggBias(dt.H0_mE[endpoint=="3 endpoints"])
ggSe.H0_mE  <- ggSe(dt.H0_mE[endpoint=="3 endpoints"])
ggCoverage.H0_mE  <- ggCoverage(dt.H0_mE[Hprojection==1])
## ggCoverage.H0_mE$data$rep

## ** export
ggsave(ggTiming.H0_mE$plot, filename = file.path(path.results,"fig-timing-H0-mE.pdf"),
       width = 10, height = 9, device = cairo_pdf)
ggsave(ggBias.H0_mE$plot, filename = file.path(path.results,"fig-bias-H0-mE.pdf"),
       device = cairo_pdf)
ggsave(ggSe.H0_mE$plot, filename = file.path(path.results,"fig-se-H0-mE.pdf"),
       device = cairo_pdf)
ggsave(ggCoverage.H0_mE$plot + coord_cartesian(ylim = c(0.9,1)),
       filename = file.path(path.results,"fig-coverage-H0-mE.pdf"),
       height = 7, width = 10, device = cairo_pdf)
saveRDS(dt.H0_mE, file = file.path(path.results,"dataCoverage.H0_mE.rds"))


## * Coverage under the alternative (multiple TTE)
cat("\n Coverage under the alternative (multiple TTE) \n")

## ** load
dt.H1_mE <- butils::sinkDirectory(path.coverageH1_mE, string.keep = "tempo")
dt.H1_mE[, endpoint := factor(endpoint,
                              levels = c("time_0.5", "timeU_0.5", "toxicity_0.5", "time_1e-12", "timeU_1e-12"),
                              labels = c("1 endpoint", "1 endpoint", "2 endpoints", "3 endpoints", "3 endpoints"))]
## dt.H1_me <- readRDS(file.path(path.results,"dataCoverage.H1_mE.rds"))
## dt.H1_mE[, endpoint := factor(endpoint,
##                               levels = c("survival (\u03C4=0.5)", "toxicity", "survival"),
##                               labels = c("1 endpoint", "2 endpoints", "3 endpoints"))]

## ** display
ggTiming.H1_mE <- ggTiming(dt.H1_mE[endpoint=="3 endpoints"], file = 1)
ggBias.H1_mE  <- ggBias(dt.H1_mE[endpoint=="3 endpoints"])
ggSe.H1_mE  <- ggSe(dt.H1_mE[endpoint=="3 endpoints"])
ggCoverage.H1_mE  <- ggCoverage(dt.H1_mE[Hprojection==1])

## ** export
ggsave(ggTiming.H1_mE$plot, filename = file.path(path.results,"fig-timing-H1-mE.pdf"),
       width = 10, height = 9, device = cairo_pdf)
ggsave(ggBias.H1_mE$plot, filename = file.path(path.results,"fig-bias-H1-mE.pdf"),
       device = cairo_pdf)
ggsave(ggSe.H1_mE$plot, filename = file.path(path.results,"fig-se-H1-mE.pdf"),
       device = cairo_pdf)
ggsave(ggCoverage.H1_mE$plot + coord_cartesian(ylim = c(0.9,1)),
       filename = file.path(path.results,"fig-coverage-H1-mE.pdf"),
       height = 7, width = 10, device = cairo_pdf)
saveRDS(dt.H1_mE, file = file.path(path.results,"dataCoverage.H1_mE.rds"))

######################################################################
### BUILD_results.R ends here
