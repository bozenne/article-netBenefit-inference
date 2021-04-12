## * Header 
## path <- "p:/Cluster/GPC/Article-inference-Ustatistic - Rao/"
## setwd(path)
## source("BATCH_SIMULATION_H1-mE-nm.R")
## sbatch -a 1-10 -J 'mytest' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_SIMULATION_H1-mE-nm.R /dev/null 

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 40}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),
                  size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * path
path <- "."
path.res <- file.path(path,"Results","SIMULATION_H1-mE-nm")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
        dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","SIMULATION_H1-mE-nm")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
        dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * libraries
library(BuyseTest)
data.table::setDTthreads(1)
library(data.table)

## * settings
n.sim <- 100
seqN <- c(30,60,120,240,400,600,1000)
BuyseTest.options(trace = 0, debug = 0)

## * job
res <- NULL
for(iN in 1:length(seqN)){ ## iN <- 5
    cat("sample size: ",seqN[iN],"\n")
    for(iSim in 1:n.sim){ ## iSim <- 1
        cat(iSim," ")
        ## ** sim data
        iN.T <- seqN[iN]*2/4
        iN.C <- seqN[iN]*2*3/4
        dt <- simBuyseTest(n.T = iN.T,
                           n.C = iN.C,
                           argsBin = list(p.T = c(0.8,0.2),
                                          p.C = c(0.2,0.8)),
                           argsTTE = list(scale.T = 2,
                                          scale.C = 1,
                                          scale.Censoring.C = 1,
                                          scale.Censoring.T = 1),
                           latent = TRUE)

        dt$One <- 1
        dt[,eventtimeCensoring:=NULL]
        setnames(dt, old = c("eventtimeUncensored","eventtime"), new = c("timeU","time"))

        ff.GS <- treatment ~ tte(timeU, status = One, threshold = 0.5) + bin(toxicity) + tte(timeU, status = One, threshold = 0)
        ff.test <- treatment ~ tte(time, status = status, threshold = 0.5) + bin(toxicity) + tte(time, status = status, threshold = 0)

        ## ** fit model
        BuyseTest.options(order.Hprojection = 2)
        tps.GS0 <- system.time(
            e.GS0 <- BuyseTest(ff.GS, data = dt, method.inference = "none")
        )
        tps.GS <- system.time(
            e.GS <- BuyseTest(ff.GS, data = dt, method.inference = "u-statistic")
        )

        tps.Gehan0 <- system.time(
            e.Gehan0 <- BuyseTest(ff.test, data = dt, scoring.rule = "Gehan", method.inference = "none")
        )
        tps.Gehan <- system.time(
            e.Gehan <- BuyseTest(ff.test, data = dt, scoring.rule = "Gehan", method.inference = "u-statistic")
        )

        BuyseTest.options(order.Hprojection = 1)
        tps.Peron0 <- system.time(
            e.Peron0 <- BuyseTest(ff.test, data = dt, scoring.rule = "Peron", method.inference = "none")
        )
        tps.Peron <- system.time(
            e.Peron <- BuyseTest(ff.test, data = dt, scoring.rule = "Peron", method.inference = "u-statistic")
        )

        ## ** confint
        iOutGS1 <- confint(e.GS, order.Hprojection = 1)
        iOutGS2 <- confint(e.GS, order.Hprojection = 2)
        iOutGehan1 <- confint(e.Gehan, order.Hprojection = 1)
        iOutGehan2 <- confint(e.Gehan, order.Hprojection = 2)
        iOutPeron1 <- confint(e.Peron)

        res <- rbind(res,
                     data.table(iOutGS1,
                                timeAll = tps.GS["elapsed"], timeEstimate = tps.GS0["elapsed"],
                                endpoint = rownames(iOutGS1), N=iN.T+iN.C, n.T = iN.T, n.C = iN.C, iter = iSim, Hprojection = 1, method = "GS"),
                     data.table(iOutGS2,
                                timeAll = tps.GS["elapsed"], timeEstimate = tps.GS0["elapsed"],
                                endpoint = rownames(iOutGS2), N=iN.T+iN.C, n.T = iN.T, n.C = iN.C, iter = iSim, Hprojection = 2, method = "GS"),
                     data.table(iOutGehan1,
                                timeAll = tps.Gehan["elapsed"], timeEstimate = tps.Gehan0["elapsed"],
                                endpoint = rownames(iOutGehan1), N=iN.T+iN.C, n.T = iN.T, n.C = iN.C, iter = iSim, Hprojection = 1, method = "Gehan"),
                     data.table(iOutGehan2,
                                timeAll = tps.Gehan["elapsed"], timeEstimate = tps.Gehan0["elapsed"],
                                endpoint = rownames(iOutGehan2), N=iN.T+iN.C, n.T = iN.T, n.C = iN.C, iter = iSim, Hprojection = 2, method = "Gehan"),
                     data.table(iOutPeron1,
                                timeAll = tps.Peron["elapsed"], timeEstimate = tps.Peron0["elapsed"],
                                endpoint = rownames(iOutPeron1), N=iN.T+iN.C, n.T = iN.T, n.C = iN.C, iter = iSim, Hprojection = 1, method = "Peron"))
    }
    cat("\n")
    saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,"(tempo).rds")))
}

## * export
saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,".rds")))

## * display
print(sessionInfo())


