### BUILD_graphtable.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 31 2021 (09:38) 
## Version: 
## Last-Updated: mar 31 2021 (10:12) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * path
## path <- "P:/Cluster/GPC/Article-inference-Ustatistic-Rao"
## setwd(path)

path.results <- "Results"
path.figures <- "figures"
path.tables <- "tables"
path.datasim <- do.call(rbind,
                        list(H0_1TTE = data.frame(folder = "SIMULATION_H0-1TTE", n = FALSE, endpoint = FALSE,
                                                  text = "balanced groups under the null with 1 TTE"), 
                             H1_1TTE = data.frame(folder = "SIMULATION_H1-1TTE", n = FALSE, endpoint = FALSE,
                                                  text = "balanced groups under the alternative with 1 TTE"), 
                             H0_mE = data.frame(folder = "SIMULATION_H0-mE", n = FALSE, endpoint = TRUE,
                                                text = "balanced groups under the null with multiple TTE"), 
                             H1_mE = data.frame(folder = "SIMULATION_H1-mE", n = FALSE, endpoint = TRUE,
                                                text = "balanced groups under the alternative with multiple TTE"), 
                             H0_1TTE_nm = data.frame(folder = "SIMULATION_H0-1TTE-nm", n = TRUE, endpoint = FALSE,
                                                     text = "unbalanced groups under the null with 1 TTE"), 
                             H1_1TTE_nm = data.frame(folder = "SIMULATION_H1-1TTE-nm", n = TRUE, endpoint = FALSE,
                                                     text = "unbalanced groups under the alternative with 1 TTE"), 
                             H0_mE_nm = data.frame(folder = "SIMULATION_H0-mE-nm", n = TRUE, endpoint = TRUE,
                                                   text = "unbalanced groups under the null with multiple TTE"),
                             H1_mE_nm = data.frame(folder = "SIMULATION_H1-mE-nm", n = TRUE, endpoint = TRUE,
                                                   text = "unbalanced groups under the alternative with multiple TTE")
                             )
                        )

## * libraries
library(data.table)
library(ggplot2)
library(ggpubr) ## neededed for graphical display
library(ggthemes) ## neededed for graphical display
library(xtable) ## neededed for display
source("FCT-gg.R")

## * Loop
n.datasim <- NROW(path.datasim)
for(iDatasim in 1:n.datasim){ ## iDatasim <- 1
    iFolder <- path.datasim[iDatasim,"folder"]

    ## ** load
    
    iData.timing <- try(readRDS(file = file.path(path.results,paste0("dataTiming-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.timing,"try-error")) next
    iData.bias <- try(readRDS(file = file.path(path.results,paste0("dataBias-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.bias,"try-error")) next
    iData.se <- try(readRDS(file = file.path(path.results,paste0("dataSe-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.se,"try-error")) next
    iData.coverage <- try(readRDS(file = file.path(path.results,paste0("dataCoverage-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.coverage,"try-error")) next
    iData.table <- try(readRDS(file = file.path(path.results,paste0("dataTable-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.table,"try-error")) next
    iData.proj <- try(readRDS(file = file.path(path.results,paste0("dataProj-",iFolder,".rds"))), silent = TRUE)
    if(inherits(iData.proj,"try-error")) next

    cat(iDatasim,") ", path.datasim[iDatasim,"text"],"\n",sep="")

    ## **  convert to figure/table
    iGGtiming <- ggTiming(iData.timing,
                          type.data = "processed", plot = FALSE) ## based on only one file i.e. a subset of the simulations

    iGGbias  <- ggBias(iData.bias,
                       type.data = "processed", plot = FALSE) ## based on only one file i.e. a subset of the simulations

    iGGse  <- ggSe(iData.se,
                   type.data = "processed", plot = FALSE) ## based on all files

    iGGcoverage  <- ggCoverage(iData.coverage,
                               type.data = "processed", plot = FALSE)

    iTable  <- createTable(iData.table,
                           digits = 3, print = FALSE, trace = FALSE, type.data = "processed",
                           by = if(path.datasim[iDatasim,"endpoint"]){"endpoint"}else{"threshold"})


    ## ** export
    iName.timing <- file.path(path.figures,paste0("figTiming-",iFolder,".pdf"))
    ggsave(iGGtiming$plot, filename = iName.timing, width = 10, height = 9, device = cairo_pdf)

    iName.bias <- file.path(path.figures,paste0("figBias-",iFolder,".pdf"))
    ggsave(iGGbias$plot, filename = iName.bias, device = cairo_pdf)

    iName.se <- file.path(path.figures,paste0("figSe-",iFolder,".pdf"))
    ggsave(iGGse$plot, filename = iName.se, device = cairo_pdf)

    iName.coverage <- file.path(path.figures,paste0("figCoverage-",iFolder,".pdf"))
    ggsave(iGGcoverage$plot + coord_cartesian(ylim = c(0.9,1)), filename = iName.coverage,
           height = 7, width = 10, device = cairo_pdf)

    fileConn <- file(file.path(path.tables,paste0("table-",iFolder,".txt")))
    writeLines(iTable$table, fileConn)
    close(fileConn)

    fileConn <- file(file.path(path.tables,paste0("tableHproj-",iFolder,".txt")))
    writeLines(capture.output(iData.proj), fileConn)
    close(fileConn)
}
######################################################################
### BUILD_graphtable.R ends here
