### BUILD_results.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 17 2018 (09:40) 
## Version: 
## Last-Updated: mar 31 2021 (09:42) 
##           By: Brice Ozenne
##     Update #: 178
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
    iPath <- file.path(path.results, iFolder)

    if(length(list.files(iPath))==0){next}
    cat(iDatasim,") ", path.datasim[iDatasim,"text"],"\n",sep="")

    ## ** load data
    iDT <- butils::sinkDirectory(iPath, string.keep = "tempo") ## , string.exclude = paste(c(2:9,0),collapse="|") 
    ## iDT <- butils::sinkDirectory(iPath, string.keep = "tempo", string.exclude = paste(c(2:9,0),collapse="|")) 

    ## ** process data
    ## n
    if(path.datasim[iDatasim,"n"]){
        iDT[, n := factor(paste0(n.T,"/",n.C), unique(paste0(n.T,"/",n.C)))]
    }
    keep.level.n <- sort(unique(iDT$n))[c(1,3,6)]
    
    ## endpoint
    if(path.datasim[iDatasim,"endpoint"]){
        iDT[, endpoint := droplevels(factor(endpoint,
                                            levels = c("time_0.5", "timeU_0.5", "toxicity_0.5", "time_1e-12", "timeU_1e-12"),
                                            labels = c("1 endpoint", "1 endpoint", "2 endpoints", "3 endpoints", "3 endpoints")))]
        iDT.bias <- iDT[endpoint=="3 endpoints"]
    }else{
        iDT[, endpoint := NULL]
        iDT.bias <- iDT[threshold==0.5]
    }

    ## ** generate data for plot and tables
    iGGtiming <- ggTiming(iDT.bias, file = 1, type.data = "raw", plot = FALSE)
    iGGbias  <- ggBias(iDT.bias, file = 1, type.data = "raw", plot = FALSE)
    iGGse  <- ggSe(iDT.bias, type.data = "raw", plot = FALSE)
    iGGcoverage  <- ggCoverage(iDT[Hprojection==1], type.data = "raw", plot = FALSE)
    iTable  <- createTable(iDT[n %in% keep.level.n & Hprojection==1], digits = 3, print = FALSE, trace = FALSE,
                           by = if(path.datasim[iDatasim,"endpoint"]){"endpoint"}else{"threshold"})

    iProj <- dcast(iDT.bias[method == "Gehan", .(rep = .N, sigma_empirical = sd(estimate),sigma_Ustat = mean(se)),
                            by = c("Hprojection","n")],
                   value.var = "sigma_Ustat",
                   formula = n + rep + sigma_empirical ~ Hprojection)

    ## ** export
    saveRDS(iGGtiming$data, file = file.path(path.results,paste0("dataTiming-",iFolder,".rds")))
    saveRDS(iGGbias$data, file = file.path(path.results,paste0("dataBias-",iFolder,".rds")))
    saveRDS(iGGse$data, file = file.path(path.results,paste0("dataSe-",iFolder,".rds")))
    saveRDS(iGGcoverage$data, file = file.path(path.results,paste0("dataCoverage-",iFolder,".rds")))
    saveRDS(iTable$data, file = file.path(path.results,paste0("dataTable-",iFolder,".rds")))
    saveRDS(iProj, file = file.path(path.results,paste0("dataProj-",iFolder,".rds")))

}


######################################################################
### BUILD_results.R ends here

