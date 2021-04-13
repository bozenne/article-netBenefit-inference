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
##    shape    scale 
## 1.289930 8.995655 
param1 <- exp(coef(AFT1))
 ##    shape     scale 
 ## 1.275269 13.765431 

AFT0.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0.cens <- exp(coef(AFT0.cens))
 ##    shape     scale 
 ## 1.369449 34.305618 
param1.cens <- exp(coef(AFT1.cens))
 ##    shape     scale 
 ## 1.490881 27.885185 

## ** toxicity: parameter of the categorical distribution
ptox.C <- prop.table(dt.prodige[,table(toxicity, bras)],2)[,"Gemcitabine"]
##           0           1           2           3           4           5 
## 0.011695906 0.029239766 0.362573099 0.391812865 0.198830409 0.005847953 
ptox.T <-prop.table(dt.prodige[,table(toxicity, bras)],2)[,"Folfirinox"]
##           0           1           2           3           4           5 
## 0.035087719 0.040935673 0.233918129 0.473684211 0.210526316 0.005847953 
## ** correlation between survival and toxicity

n <- 1e4
seqRho <- seq(-0.1,0,by=0.01)

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
##       rho    delta1      delta2     Delta         sd1          sd2          sd
##  1: -0.10 0.2574496 -0.03998591 0.2174637 0.002619748 0.0004570351 0.002678611
##  2: -0.09 0.2626466 -0.03729770 0.2253489 0.001217537 0.0003174556 0.001166779
##  3: -0.08 0.2614341 -0.03608493 0.2253491 0.001510400 0.0003817177 0.001691942
##  4: -0.07 0.2619730 -0.03350018 0.2284729 0.002100621 0.0003566316 0.002224376
##  5: -0.06 0.2571747 -0.03122962 0.2259450 0.002492393 0.0004493182 0.002686840
##  6: -0.05 0.2615422 -0.02826779 0.2332744 0.001441908 0.0003612850 0.001199878
##  7: -0.04 0.2628981 -0.02557217 0.2373259 0.002401048 0.0004289239 0.002369637
##  8: -0.03 0.2611181 -0.02140634 0.2397118 0.003070295 0.0005792831 0.002752112
##  9: -0.02 0.2600821 -0.01769573 0.2423864 0.003378429 0.0003689707 0.003464852
## 10: -0.01 0.2593785 -0.01425307 0.2451254 0.002013681 0.0005714573 0.002134403
## 11:  0.00 0.2563715 -0.01078417 0.2455874 0.002844241 0.0005401532 0.003192720

## ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
## ##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## ## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1583;0.3982] 1.3650e-05 ***

Mdiff <- cbind(rho = dtS.out$rho,sweep(dtS.out[,c("delta1","delta2","Delta")], FUN = "-", STATS = delta.obs, MARGIN = 2))
cbind(Mdiff, absError = abs(Mdiff$delta1)+abs(Mdiff$delta2))
##      rho      delta1        delta2       Delta   absError
## 1  -0.10 -0.04643492 -0.0187635416 -0.06519847 0.06519847
## 2  -0.09 -0.04123788 -0.0160753336 -0.05731321 0.05731321
## 3  -0.08 -0.04245043 -0.0148625666 -0.05731300 0.05731300
## 4  -0.07 -0.04191147 -0.0122778106 -0.05418928 0.05418928
## 5  -0.06 -0.04670984 -0.0100072526 -0.05671709 0.05671709
## 6  -0.05 -0.04234227 -0.0070454206 -0.04938769 0.04938769
## 7  -0.04 -0.04098641 -0.0043497996 -0.04533621 0.04533621
## 8  -0.03 -0.04276642 -0.0001839686 -0.04295039 0.04295039
## 9  -0.02 -0.04380237  0.0035266334 -0.04027574 0.04732900
## 10 -0.01 -0.04450600  0.0069692984 -0.03753670 0.05147530
## 11  0.00 -0.04751299  0.0104382014 -0.03707479 0.05795119

## * Power calculation
library(BuyseTest)
library(data.table)

simTrial <- function(n.C,n.T){
    d <- simBuyseTest(n.C = n.C,
                      n.T = n.T,
                      argsBin = list(p.C = c(0.01170, 0.02924, 0.36257, 0.39181, 0.19883, 0.00585 ), ## round(ptox.C,5)
                                     p.T = c(0.03508, 0.04094, 0.23392, 0.47368, 0.21053, 0.00585 ), ## round(ptox.T,5)
                                     rho.T = -0.03, rho.C = -0.03),
                      argsTTE = list(scale.C = 8.995655 , ## param0["scale"]
                                     scale.T = 13.76543 , ## param1["scale"]
                                     shape.C = 1.28993, ## param0["shape"]
                                     shape.T = 1.275269, ## param1["shape"]
                                     scale.Censoring.C = 34.30562,  ## param0.cens["scale"]
                                     scale.Censoring.T = 27.88519, ## param1.cens["scale"]
                                     shape.Censoring.C = 1.369449, ## param0.cens["shape"]
                                     shape.Censoring.T = 1.490881) ## param1.cens["shape"]
                      )
    d$toxicity.num <- as.numeric(d$toxicity)
    return(d)
}
## BuyseTest(treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"), data = simTrial(1000,1000))


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

 ## - statistic   : net benefit (null hypothesis Delta=0)
 ##     endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
 ##    eventtime         2  50  50         0.259      0.1105  0.1119          0.602
 ##                        100 100        0.2606      0.0785  0.0791           0.89
 ##                        150 150        0.2603      0.0631  0.0646          0.975
 ##                        200 200        0.2602      0.0532  0.0559          0.995
 ## toxicity.num       0.5  50  50        0.2396      0.1128  0.1149          0.522
 ##                        100 100        0.2407      0.0804  0.0812           0.83
 ##                        150 150        0.2402      0.0647  0.0663          0.942
 ##                        200 200        0.2403      0.0547  0.0574          0.984

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

tps0
##   user  system elapsed 
## 22.045   0.004  22.058 

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

 ## - statistic   : net benefit (null hypothesis Delta=0)
 ##     endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
 ##    eventtime         2  50  50         0.259      0.1105  0.1119          0.602
 ##                        100 100        0.2606      0.0785  0.0791           0.89
 ##                        150 150        0.2603      0.0631  0.0646          0.975
 ##                        200 200        0.2602      0.0532  0.0559          0.995
 ## toxicity.num       0.5  50  50        0.2396      0.1128  0.1149          0.522
 ##                        100 100        0.2407      0.0804  0.0812           0.83
 ##                        150 150        0.2402      0.0647  0.0663          0.942
 ##                        200 200        0.2403      0.0547  0.0574          0.984

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
                               sample.size = c(50,92,93,94,95,96,100),
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

