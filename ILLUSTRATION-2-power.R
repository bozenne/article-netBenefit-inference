## * data + packages
source("ILLUSTRATION-0-data-management.R")

## * Analysis
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 0.5, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic", trace = FALSE)
delta.obs <- c(delta1 = as.double(coef(e.BT)[1]),
               delta2 = as.double(diff(coef(e.BT))),
               Delta = as.double(coef(e.BT)[2]))
 ##     delta1      delta2       Delta 
 ## 0.30388451 -0.02122237  0.28266214 

## * Estimation of the simulation parameters

## ** survival: parameters of the weibull distribution
AFT0 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0 <- exp(coef(AFT0))
param1 <- exp(coef(AFT1))

AFT0.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0.cens <- exp(coef(AFT0.cens))
param1.cens <- exp(coef(AFT1.cens))

## ** toxicity: parameter of the categorical distribution
ptox.C <- c("0" = 1.2/100,
            "1" = 2.9/100,
            "2" = 36.2/100,
            "3" = 39.2/100,
            "4" = 19.9/100,
            "5" = 0.6/100)
ptox.T <-  c("0" = 3.5/100,
             "1" = 4.1/100,
             "2" = 23.4/100,
             "3" = 47.3/100,
             "4" = 21.1/100,
             "5" = 0.6/100)

## ** correlation between survival and toxicity
n <- 5e3
seqRho <- seq(-0.3,0.3,by=0.025)

warperCor <- function(iProc){
    argsTTE <- list(scale.C = param0["scale"],
                    scale.T = param1["scale"],
                    shape.C = param0["shape"],
                    shape.T = param1["shape"],
                    scale.Censoring.C = 1e5*param0.cens["scale"], ## no censoring
                    scale.Censoring.T = 1e5*param1.cens["scale"], ## no censoring
                    shape.Censoring.C = param0.cens["shape"],
                    shape.Censoring.T = param1.cens["shape"])

    out <- NULL
    for(iRho in 1:length(seqRho)){ ## iRho <- 1
        argsBin <- list(p.C = ptox.C,
                        p.T = ptox.T,
                        rho.T = seqRho[iRho], rho.C = seqRho[iRho])


        ## simulate data for a given correlation
        dt.sim <- simBuyseTest(n,
                               argsBin = argsBin,
                               argsTTE = argsTTE,
                               latent = TRUE)
        dt.sim$toxicity.num <- as.numeric(dt.sim$toxicity)-1

        ## estimate net benefit
        iE.BT <- BuyseTest(treatment ~ tte(eventtimeUncensored, status = status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"),
                           data = dt.sim, trace = FALSE)
        ## summary(iE.BT)
        ## iE.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, threshold = 0.5), data = dt.sim, trace = FALSE)
        ## summary(iE.BT2)
        out <- rbind(out,
                     cbind(rep = iProc, rho = seqRho[iRho], delta1 = coef(iE.BT)[1], delta2 = coef(iE.BT)[2]-coef(iE.BT)[1], Delta = coef(iE.BT)[2])
                     )
    }
    return(out)
}

set.seed(10)
cpus <- 10
n.sim <- 10

cl <- snow::makeSOCKcluster(cpus)
doSNOW::registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style=3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

parallel::clusterExport(cl, varlist = c("param0","param1","param0.cens","param1.cens","seqRho"))

ls.res <- foreach::`%dopar%`(
                       foreach::foreach(i=1:n.sim, .options.snow=opts), {
                           require(BuyseTest)
                           warperCor(i)
                       })
parallel::stopCluster(cl)

dt.out <- as.data.table(do.call(rbind,ls.res))
## saveRDS(dt.out, file = file.path("Results","rho-SIMULATION-power.rds"))

dtS.out <- dt.out[, .(delta1 = mean(delta1), delta2 = mean(delta2), Delta = mean(Delta),
                      sd1 = sd(delta1)/sqrt(.N), sd2 = sd(delta2)/sqrt(.N), sd = sd(Delta)/sqrt(.N)), by = "rho"]
dtS.out
##               rho    delta1       delta2     Delta         sd1          sd2          sd
##  1: -3.000000e-01 0.2589394 -0.061794032 0.1971454 0.002273165 0.0003896854 0.002415783
##  2: -2.750000e-01 0.2642665 -0.060482268 0.2037843 0.003930216 0.0005924597 0.004162672
##  3: -2.500000e-01 0.2595744 -0.058676296 0.2008981 0.003549636 0.0002764130 0.003556239
##  4: -2.250000e-01 0.2650715 -0.055129272 0.2099422 0.003952291 0.0005481797 0.004018452
##  5: -2.000000e-01 0.2582152 -0.054083252 0.2041320 0.004356295 0.0008357436 0.004583327
##  6: -1.750000e-01 0.2543116 -0.051381948 0.2029296 0.004279700 0.0007828054 0.004004096
##  7: -1.500000e-01 0.2574692 -0.048533316 0.2089359 0.002980460 0.0004533307 0.003099940
##  8: -1.250000e-01 0.2576633 -0.045158592 0.2125047 0.002613338 0.0004717178 0.002756413
##  9: -1.000000e-01 0.2569524 -0.039703896 0.2172485 0.004550971 0.0007612167 0.004452540
## 10: -7.500000e-02 0.2638728 -0.034684316 0.2291885 0.003018895 0.0008868919 0.002934202
## 11: -5.000000e-02 0.2599286 -0.028110652 0.2318180 0.002846804 0.0005439338 0.002629213
## 12: -2.500000e-02 0.2586773 -0.019501576 0.2391757 0.003901417 0.0008683635 0.004345399
## 13:  5.551115e-17 0.2577217 -0.010588204 0.2471335 0.002897286 0.0004456269 0.002934284
## 14:  2.500000e-02 0.2634448 -0.001434300 0.2620105 0.002469265 0.0008201922 0.002176349
## 15:  5.000000e-02 0.2630703  0.005811684 0.2688820 0.002397104 0.0006247747 0.002202936
## 16:  7.500000e-02 0.2558377  0.012298328 0.2681360 0.002788943 0.0007143914 0.002798140
## 17:  1.000000e-01 0.2622366  0.016729088 0.2789657 0.002003683 0.0005780145 0.001955335
## 18:  1.250000e-01 0.2629381  0.020908088 0.2838461 0.001923403 0.0004420016 0.001867996
## 19:  1.500000e-01 0.2547908  0.023580344 0.2783711 0.001899027 0.0006119372 0.001506621
## 20:  1.750000e-01 0.2639497  0.026431024 0.2903808 0.003785003 0.0006001119 0.003400504
## 21:  2.000000e-01 0.2569959  0.028818584 0.2858145 0.003677529 0.0006210232 0.003441335
## 22:  2.250000e-01 0.2611091  0.030769944 0.2918791 0.002882215 0.0005496446 0.002659747
## 23:  2.500000e-01 0.2670079  0.032936464 0.2999444 0.002883019 0.0006817218 0.002400987
## 24:  2.750000e-01 0.2619929  0.033633416 0.2956263 0.004140570 0.0005118157 0.003876921
## 25:  3.000000e-01 0.2648798  0.034228744 0.2991085 0.002373132 0.0004103843 0.002522153
##               rho    delta1       delta2     Delta         sd1          sd2          sd

## ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
## ##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## ## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1583;0.3982] 1.3650e-05 ***

Mdiff <- cbind(rho = dtS.out$rho,sweep(dtS.out[,c("delta1","delta2","Delta")], FUN = "-", STATS = delta.obs, MARGIN = 2))
cbind(Mdiff, absError = abs(Mdiff$delta1)+abs(Mdiff$delta2))
##              rho      delta1       delta2        Delta   absError
## 1  -3.000000e-01 -0.04494507 -0.040571665 -0.085516735 0.08551673
## 2  -2.750000e-01 -0.03961797 -0.039259901 -0.078877867 0.07887787
## 3  -2.500000e-01 -0.04431009 -0.037453929 -0.081764023 0.08176402
## 4  -2.250000e-01 -0.03881301 -0.033906905 -0.072719915 0.07271991
## 5  -2.000000e-01 -0.04566926 -0.032860885 -0.078530143 0.07853014
## 6  -1.750000e-01 -0.04957293 -0.030159581 -0.079732511 0.07973251
## 7  -1.500000e-01 -0.04641527 -0.027310949 -0.073726219 0.07372622
## 8  -1.250000e-01 -0.04622122 -0.023936225 -0.070157443 0.07015744
## 9  -1.000000e-01 -0.04693206 -0.018481529 -0.065413591 0.06541359
## 10 -7.500000e-02 -0.04001173 -0.013461949 -0.053473683 0.05347368
## 11 -5.000000e-02 -0.04395587 -0.006888285 -0.050844155 0.05084415
## 12 -2.500000e-02 -0.04520718  0.001720791 -0.043486391 0.04692797
## 13  5.551115e-17 -0.04616278  0.010634163 -0.035528615 0.05679694
## 14  2.500000e-02 -0.04043969  0.019788067 -0.020651623 0.06022776
## 15  5.000000e-02 -0.04081421  0.027034051 -0.013780163 0.06784827
## 16  7.500000e-02 -0.04804679  0.033520695 -0.014526095 0.08156749
## 17  1.000000e-01 -0.04164791  0.037951455 -0.003696459 0.07959937
## 18  1.250000e-01 -0.04094645  0.042130455  0.001184009 0.08307690
## 19  1.500000e-01 -0.04909374  0.044802711 -0.004291027 0.09389645
## 20  1.750000e-01 -0.03993477  0.047653391  0.007718617 0.08758817
## 21  2.000000e-01 -0.04688857  0.050040951  0.003152381 0.09692952
## 22  2.250000e-01 -0.04277540  0.051992311  0.009216913 0.09476771
## 23  2.500000e-01 -0.03687657  0.054158831  0.017282257 0.09103541
## 24  2.750000e-01 -0.04189161  0.054855783  0.012964177 0.09674739
## 25  3.000000e-01 -0.03900474  0.055451111  0.016446369 0.09445585


## * Power calculation
library(BuyseTest)
library(data.table)


simTrial <- function(n.C,n.T){
    d <- simBuyseTest(n.C = n.C,
                      n.T = n.T,
                      argsBin = list(p.C = c("0" = 1.2/100,
                                             "1" = 2.9/100,
                                             "2" = 36.2/100,
                                             "3" = 39.2/100,
                                             "4" = 19.9/100,
                                             "5" = 0.6/100),
                                     p.T = c("0" = 3.5/100,
                                             "1" = 4.1/100,
                                             "2" = 23.4/100,
                                             "3" = 47.3/100,
                                             "4" = 21.1/100,
                                             "5" = 0.6/100),
                                     rho.T = -0.025, rho.C = -0.025),
                      argsTTE = list(scale.C = 8.55,
                                     scale.T = 13.35,
                                     shape.C = 1.29,
                                     shape.T = 1.28,
                                     scale.Censoring.C = 34.31,
                                     scale.Censoring.T = 27.88,
                                     shape.Censoring.C = 1.37,
                                     shape.Censoring.T = 1.49)
                      )
    d$toxicity.num <- as.numeric(d$toxicity)
    return(d)
}
## BuyseTest(treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"), data = simTrial(500,500))


## ** 1 CPUS, n.rep = 100
tps0 <- system.time(
    e0.power <- powerBuyseTest(sim = simTrial,
                               formula = treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"),
                               sample.size = c(50,100,150,200),
                               n.rep = 100, seed = 10, 
                               method.inference = "u-statistic")
)

summary(e0.power, endpoint = c("eventtime_2","toxicity.num_0.5"))
##         Simulation study with Generalized pairwise comparison
##         with 100 samples

##  - statistic   : net benefit (null hypothesis Delta=0)
##      endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##     eventtime         2  50  50        0.2797      0.1247  0.1105           0.65
##                         100 100        0.2676      0.0876  0.0785           0.88
##                         150 150        0.2663      0.0669  0.0641           0.99
##                         200 200        0.2657      0.0568  0.0556              1
##  toxicity.num       0.5  50  50        0.2616      0.1254  0.1138            0.6
##                         100 100        0.2499      0.0878  0.0806           0.81
##                         150 150        0.2485      0.0672  0.0658           0.98
##                         200 200        0.2477      0.0574  0.0571           0.99

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

tps0
 ##   user  system elapsed 
 ## 23.247   0.008  23.257 

## ** 10 CPUS, n.rep = 1000
tps1 <- system.time(
    e1.power <- powerBuyseTest(sim = simTrial,
                               formula = treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"),
                               sample.size = c(50,100,150,200),
                               n.rep = 1000, seed = 10, cpus = 10,
                               method.inference = "u-statistic")
)

summary(e1.power, endpoint = c("eventtime_2","toxicity.num_0.5"))
##         Simulation study with Generalized pairwise comparison
##         with 1000 samples

##  - statistic   : net benefit (null hypothesis Delta=0)
##      endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##     eventtime         2  50  50          0.27      0.1108  0.1109          0.642
##                         100 100        0.2718       0.078  0.0784          0.917
##                         150 150        0.2723      0.0623   0.064          0.985
##                         200 200        0.2721      0.0531  0.0554          0.997
##  toxicity.num       0.5  50  50        0.2507      0.1136  0.1141          0.551
##                         100 100        0.2525        0.08  0.0806           0.87
##                         150 150         0.253       0.064  0.0658          0.964
##                         200 200         0.253      0.0547   0.057          0.992

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

tps1
##  user  system elapsed 
## 0.588   0.064  44.410 

## ** refinement
tps2 <- system.time(
    e2.power <- powerBuyseTest(sim = simTrial,
                               formula = treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"),
                               sample.size = c(50,83,84,85,86,87,100),
                               n.rep = 1000, seed = 10, cpus = 10,
                               method.inference = "u-statistic")
)

summary(e2.power)
##         Simulation study with Generalized pairwise comparison
##         with 1000 samples

##  - statistic   : net benefit (null hypothesis Delta=0)
##      endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##  toxicity.num       0.5  50  50        0.2516      0.1161   0.114          0.571
##                          83  83         0.254      0.0893  0.0884          0.791
##                          84  84        0.2541      0.0888  0.0879            0.8
##                          85  85        0.2541      0.0887  0.0874          0.803
##                          86  86        0.2542      0.0882  0.0868            0.8
##                          87  87        0.2538      0.0877  0.0864          0.804
##                         100 100        0.2538      0.0813  0.0805          0.864

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

