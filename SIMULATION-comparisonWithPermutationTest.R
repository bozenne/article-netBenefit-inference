### ILLUSTRATION-simulations.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  3 2020 (13:22) 
## Version: 
## Last-Updated: okt  8 2020 (17:31) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(BuyseTest)
library(ggplot2)
library(pbapply)

## * Simulation study: Gehan under the null
## ** Type 1 error: permutation
warperPerm <- function(i, m, n){ ## total sample size 8*n
    mu <- 0.3
    sigma <- c(0.06,0.3)
    df <- rbind(data.frame(estradiol = rnorm(2*m, mean = mu, sd = sigma[1]),
                           limit = c(rep(0.04, m), rep(0.09, m)),
                           group = "C"),
                data.frame(estradiol = rnorm(2*n, mean = mu, sd = sigma[2]),
                           limit = c(rep(0.04, n), rep(0.09, n)),
                           group = "T"))
    df$estradiol.obs <- pmax(df$estradiol, df$limit)
    df$observed <- as.numeric(df$estradiol >= df$limit)
    
    iBT <- BuyseTest(group ~ tte(estradiol.obs, status = observed, threshold = 0.05, censoring = "left"), 
                     data = df, scoring.rule = "Gehan", trace = 0,
                     method.inference = "studentized permutation", n.resampling = 1e3)

    iBT2 <- iBT
    iBT2@method.inference[] <- "u-statistic"
    attr(iBT2@method.inference,"permutation") <- FALSE
    attr(iBT2@method.inference,"ustatistic") <- TRUE
    
    ## Perm
    iOut <- rbind(data.frame(i = i, inference = "permutation", confint(iBT, method.ci.resampling  = "percentile")),
                  data.frame(i = i, inference = "studentized permutation", confint(iBT, method.ci.resampling  = "studentized")),
                  data.frame(i = i, inference = "u-statistic", confint(iBT2)))
    rownames(iOut) <- NULL
    return(iOut)
}

warperPerm(i = 1, m = 100, n = 100)

cpus <- 32
cl <- snow::makeSOCKcluster(cpus)
doSNOW::registerDoSNOW(cl)

set.seed(10)
n.sim <- 10000
pb <- txtProgressBar(max = n.sim, style=3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
ls.res <- foreach::`%dopar%`(
                       foreach::foreach(i=1:n.sim, .options.snow=opts, .packages = c("data.table","BuyseTest")), {
                           ## iRes <- NULL
                           ## or(iN in c(30,50,100)){
                           ##     iRes <- rbind(iRes,
                           ##                   cbind(i=i,n=iN,warperPerm(i=i,m=iN,n=iN)))
                           ## }
                           iRes <- rbind(cbind(i=i,N=30+30,warperPerm(i=i,m=30,n=30)),
                                         cbind(i=i,N=11+33,warperPerm(i=i,m=11,n=33)))
                           return(iRes)
                       })

dt.res <- as.data.table(do.call(rbind,ls.res))
dt.res[,.(rep = .N, estimate = mean(estimate), type1error = mean(p.value<=0.05)), by = c("inference","N")]
##                  inference  N   rep      estimate type1error
## 1:             permutation 60 10000  4.202778e-05     0.0865
## 2: studentized permutation 60 10000  4.202778e-05     0.0498
## 3:             u-statistic 60 10000  4.202778e-05     0.0478
## 4:             permutation 44 10000 -4.777548e-04     0.0154
## 5: studentized permutation 44 10000 -4.777548e-04     0.0476
## 6:             u-statistic 44 10000 -4.777548e-04     0.0497

## * [NO USED] Bias: t-test vs.GPC
warper <- function(i, n = c(4,7), m = c(12,21), return.data = FALSE){
                                        # ## n <- c(4,7)*500; m <- c(12,21)*500;
    mu <- 0.3
    sigma <- c(0.06,0.3)
    df <- rbind(data.frame(estradiol = rnorm(sum(m), mean = mu, sd = sigma[1]),
                           limit = c(rep(0.04, m[1]), rep(0.09, m[2])),
                           group = "C"),
                data.frame(estradiol = rnorm(sum(n), mean = mu, sd = sigma[2]),
                           limit = c(rep(0.04, n[1]), rep(0.09, n[2])),
                           group = "T"))
    df$estradiol.obs <- pmax(df$estradiol, df$limit)
    df$observed <- as.numeric(df$estradiol >= df$limit)
    ## as.data.table(df)[,mean(observed), by = "group"]

    ## ggplot(df, aes(estradiol, fill = group)) + geom_histogram(position="identity", alpha = 0.7)
    ## ggplot(df, aes(estradiol.obs, fill = group)) + geom_histogram(position="identity", alpha = 0.7)
    
    iBT <- BuyseTest(group ~ tte(estradiol.obs, status = observed, threshold = 0.05, censoring = "left"), 
                     data = df, scoring.rule = "Gehan", trace = 0)
    tau <- iBT@threshold
    Uplus <- pnorm(-tau/sqrt(sigma[1]^2+sigma[2]^2))
    iTT <- t.test(estradiol.obs ~ group, df)
    iOut <- rbind(data.frame(i = i, statistic = "netBenefit", confint(iBT, statistic = "netBenefit", null = 0)),
                  data.frame(i = i, statistic = "favorable", confint(iBT, statistic = "favorable", null = Uplus)),
                  data.frame(i = i, statistic = "unfavorable", confint(iBT, statistic = "unfavorable", null = Uplus)),
                  data.frame(i = i, statistic = "t-test", estimate = diff(iTT$estimate), se = NA, lower.ci = iTT$conf.int[1], upper.ci = iTT$conf.int[2], p.value = iTT$p.value))
    rownames(iOut) <- NULL
    if(return.data){
        attr(iOut,"data") <- as.data.table(df)
    }
    return(iOut)
}

## run simulation
set.seed(10)
dt <- as.data.table(do.call(rbind,pblapply(1:2500, warper)))

## run large sample values
set.seed(10)
out <- warper(1, n = c(500,500), m = c(1500, 1500), return.data = TRUE)
out
## :   i   statistic    estimate         se    lower.ci    upper.ci      p.value
## : 1 1  netBenefit -0.04069532 0.02664972 -0.09276942  0.01160081 1.271697e-01
## : 2 1   favorable  0.41722424 0.01364221  0.39075362  0.44418070 1.928314e-01
## : 3 1 unfavorable  0.45791956 0.01360429  0.43140101  0.48467870 9.182646e-02
## : 4 1      t-test  0.03161912         NA -0.04630887 -0.01692937 2.599459e-05

## > data.gg[,mean(observed),by="group"]
##    group        V1
## 1:     C 0.9996667
## 2:     T 0.7790000

## display distributions
data.gg <- melt(attr(out,"data"), id.vars = c("limit","group","observed"),
                value.name = "estradiol", variable.name = "origin")
data.gg[, origin := factor(origin, levels = c("estradiol","estradiol.obs"), labels = c("simulated","observed"))]
data.gg[, group := factor(group, levels = c("C","T"), labels = c("non-user","user"))]

gg <- ggplot(data.gg, aes(estradiol, fill = group))
gg <- gg + geom_histogram(aes(y = stat(count) / sum(count)), position="identity", alpha = 0.7)
gg <- gg + facet_grid(~origin)  + ylab("relative frequency")



## * [NOT USED] Inflated type 1 error of the permutation test (example from the litterature)
n <- 1e2
warperPermLitt <- function(i, n){
    df <- rbind(data.frame(value = rnorm(n, mean = 0, sd = 5), group = "C"),
                data.frame(value = rnorm(n, mean = 0, sd = 1), group = "T")
                )
    e1.BT <- BuyseTest(group ~ cont(value), data = df,
                       method.inference = "studentized permutation", trace = 0)
    e2.BT <- e1.BT
    e2.BT@method.inference[] <- "u-statistic"
    attr(e2.BT@method.inference,"permutation") <- FALSE
    attr(e2.BT@method.inference,"ustatistic") <- TRUE
    
    return(rbind(data.table(i = i, inference = "permutation", confint(e1.BT, method.ci.resampling = "percentile")),
                 data.table(i = i, inference = "studentized permutation", confint(e1.BT, method.ci.resampling = "studentized")),
                 data.table(i = i, inference = "u-statistic", confint(e2.BT)))
           )
}
cpus <- 4
cl <- snow::makeSOCKcluster(cpus)
doSNOW::registerDoSNOW(cl)

set.seed(10)
n.sim <- 5e3
pb <- txtProgressBar(max = n.sim, style=3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
ls.res <- foreach::`%dopar%`(
    foreach::foreach(i=1:n.sim, .options.snow=opts, .packages = c("data.table","BuyseTest")), {
        warperPermLitt(i = i, n = n)
    })

dt.res <- as.data.table(do.call(rbind,ls.res))
dt.res[,.(type1error = mean(p.value<=0.05)), by = "method"]

## :                     method type1error
## : 1:             permutation     0.0784
## : 2: studentized permutation     0.0448
## : 3:             u-statistic     0.0450

##----------------------------------------------------------------------
### ILLUSTRATION-simulations.R ends here
