## * data + packages
source("ILLUSTRATION-data-management.R")

## * Analysis
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 0.5, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic", trace = FALSE)
delta.obs <- c(delta1 = as.double(coef(e.BT)[1]),
               delta2 = as.double(diff(coef(e.BT))),
               Delta = as.double(coef(e.BT)[2]))

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

set.seed(10)
out <- NULL
for(iRho in 1:length(seqRho)){ ## iRho <- 1
    cat(iRho, " ")
    argsBin <- list(p.C = ptox.C,
                    p.T = ptox.T,
                    rho.T = seqRho[iRho], rho.C = seqRho[iRho])

    argsTTE <- list(scale.C = param0["scale"],
                    scale.T = param1["scale"],
                    shape.C = param0["shape"],
                    shape.T = param1["shape"],
                    scale.Censoring.C = param0.cens["scale"],
                    scale.Censoring.T = param1.cens["scale"],
                    shape.Censoring.C = param0.cens["shape"],
                    shape.Censoring.T = param1.cens["shape"])

    for(iRep in 1:10){
    ## simulate data for a given correlation
    dt.sim <- simBuyseTest(n,
                           argsBin = argsBin,
                           argsTTE = argsTTE,
                           latent = TRUE)
    dt.sim$toxicity.num <- as.numeric(dt.sim$toxicity)-1
    dt.sim$one <- 1

        ## estimate net benefit
        iE.BT <- BuyseTest(treatment ~ tte(eventtimeUncensored, status = one, threshold = 2) + cont(toxicity.num, threshold = 0.5, operator = "<0"),
                           data = dt.sim, trace = FALSE)
        ## summary(iE.BT)
        ## iE.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, threshold = 0.5), data = dt.sim, trace = FALSE)
        ## summary(iE.BT2)
        out <- rbind(out,
                     cbind(rep = iRep, rho = seqRho[iRho], delta1 = coef(iE.BT)[1], delta2 = coef(iE.BT)[2]-coef(iE.BT)[1], Delta = coef(iE.BT)[2])
                     )
    }
}
rownames(out) <- NULL
dt.out <- as.data.table(out)
## saveRDS(dt.out, file = file.path("results","power-select-rho.rds"))

dtS.out <- dt.out[, .(delta1 = mean(delta1), delta2 = mean(delta2), Delta = mean(Delta),
                      sd1 = sd(delta1)/sqrt(.N), sd2 = sd(delta2)/sqrt(.N), sd = sd(Delta)/sqrt(.N)), by = "rho"]
dtS.out
##               rho    delta1       delta2     Delta         sd1          sd2          sd
##  1: -3.000000e-01 0.3021378 -0.070007184 0.2321306 0.003498308 0.0004319403 0.003708506
##  2: -2.750000e-01 0.3118252 -0.068629264 0.2431959 0.003664251 0.0004966306 0.003841853
##  3: -2.500000e-01 0.3064834 -0.066350784 0.2401326 0.002411327 0.0005893630 0.002670403
##  4: -2.250000e-01 0.3004858 -0.064570656 0.2359151 0.004659301 0.0006974643 0.004606008
##  5: -2.000000e-01 0.3046391 -0.061717904 0.2429212 0.002540626 0.0004493292 0.002866688
##  6: -1.750000e-01 0.3060230 -0.058884680 0.2471383 0.003449632 0.0003792834 0.003394604
##  7: -1.500000e-01 0.3050660 -0.054915276 0.2501507 0.003222949 0.0004324646 0.002997364
##  8: -1.250000e-01 0.2993711 -0.050421732 0.2489494 0.003816327 0.0007258500 0.004259065
##  9: -1.000000e-01 0.3045447 -0.045628644 0.2589160 0.002838737 0.0006633343 0.002710679
## 10: -7.500000e-02 0.3051495 -0.038965756 0.2661838 0.003796162 0.0005716145 0.003647905
## 11: -5.000000e-02 0.3064056 -0.030019148 0.2763864 0.003812348 0.0006218950 0.004300840
## 12: -2.500000e-02 0.3045620 -0.020427152 0.2841348 0.002106051 0.0007286691 0.002450817
## 13:  5.551115e-17 0.3039484 -0.010659288 0.2932891 0.002855730 0.0006988743 0.002784909
## 14:  2.500000e-02 0.3063017 -0.000640940 0.3056607 0.003035527 0.0006574316 0.003563369
## 15:  5.000000e-02 0.3050298  0.008096704 0.3131265 0.003975522 0.0007395209 0.004079234
## 16:  7.500000e-02 0.3054514  0.016046928 0.3214983 0.001494886 0.0006442834 0.001568485
## 17:  1.000000e-01 0.2997560  0.020995112 0.3207511 0.004105756 0.0008039558 0.003901799
## 18:  1.250000e-01 0.3001843  0.025575488 0.3257598 0.003502773 0.0004666966 0.003368891
## 19:  1.500000e-01 0.3011582  0.029915788 0.3310740 0.004943213 0.0004720071 0.005065082
## 20:  1.750000e-01 0.2994022  0.033510880 0.3329131 0.004554030 0.0005385758 0.004377399
## 21:  2.000000e-01 0.2979646  0.036381696 0.3343463 0.002934777 0.0006603242 0.002767582
## 22:  2.250000e-01 0.3024766  0.038860136 0.3413367 0.003526384 0.0009864679 0.003816511
## 23:  2.500000e-01 0.3063496  0.039675780 0.3460254 0.003403475 0.0003694940 0.003242897
## 24:  2.750000e-01 0.3083227  0.041779072 0.3501018 0.003328045 0.0008551262 0.002956948
## 25:  3.000000e-01 0.3011541  0.044098280 0.3452524 0.001860197 0.0003915302 0.001836559
##               rho    delta1       delta2     Delta         sd1          sd2          sd
sweep(dtS.out[,c("delta1","delta2","Delta")], FUN = "-", STATS = delta.obs, MARGIN = 2)


## ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
## ##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## ## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1583;0.3982] 1.3650e-05 ***

## reasonnable delta2 are in the range [0.09-0.21] (tolerance 0.01)
## all delta1 are too low but the highest in the range [0.09-0.21] are obtained for 0.09 and 0.18
## 0.18 is kept but this is a bit arbitrary

## * Power calculation
library(BuyseTest)
library(data.table)

system.time(
    e.power <- powerBuyseTest(sim = function(n.C,n.T){
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
                                         rho.T = 0.18, rho.C = 0.18),
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
    },
    formula = treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5),
    sample.size = c(50,100,150,200),
    n.rep = 1000, seed = 10, cpus = 4,
    method.inference = "u-statistic")
)n


## ** 1 CPUS, n.rep = 100
summary(e.power, endpoint = c("eventtime_2","toxicity.num_0.5"))
## > summary(e.power, endpoint = c("eventtime_2","toxicity.num_0.5"))
##         Simulation study with Generalized pairwise comparison
##         with 100 samples

##  > statistic   : net benefit (null hypothesis Delta=0)
##      endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##     eventtime         2  50  50        0.2797      0.1247  0.1105           0.65
##                         100 100        0.2676      0.0876  0.0785           0.88
##                         150 150        0.2663      0.0669  0.0641           0.99
##                         200 200        0.2657      0.0568  0.0556              1
##  toxicity.num       0.5  50  50         0.254       0.131  0.1139           0.57
##                         100 100        0.2424      0.0914  0.0808           0.79
##                         150 150        0.2407      0.0699   0.066           0.94
##                         200 200        0.2398      0.0586  0.0573           0.99

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

 ##   user  system elapsed 
 ## 125.37    2.22  131.34 

## ** 4 CPUS, n.rep = 1000
summary(e.power, endpoint = c("eventtime_2","toxicity.num_0.5"))
   ## user  system elapsed 
   ## 1.49    0.21  263.48 

## Simulation study with Generalized pairwise comparison
## with 1000 samples

##  > statistic   : net benefit (null hypothesis Delta=0)
##      endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
##     eventtime         2  50  50        0.2709      0.1113  0.1108          0.659
##                         100 100        0.2757      0.0762  0.0783          0.929
##                         150 150        0.2757      0.0627  0.0639           0.99
##                         200 200        0.2744      0.0553  0.0554          0.998
##  toxicity.num       0.5  50  50        0.2455      0.1138  0.1144          0.564
##                         100 100        0.2509      0.0783  0.0807          0.866
##                         150 150        0.2508      0.0644  0.0659           0.96
##                         200 200        0.2497      0.0568  0.0571          0.993

##  n.T          : number of observations in the treatment group
##  n.C          : number of observations in the control group
##  mean.estimate: average estimate over simulations
##  sd.estimate  : standard deviation of the estimate over simulations
##  mean.se      : average estimated standard error of the estimate over simulations
##  rejection    : frequency of the rejection of the null hypothesis over simulations
## (standard error: H-projection of order 1| p-value: after transformation) 

## ** refinement
system.time(
    e.power <- powerBuyseTest(sim = function(n.C,n.T){
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
                                         rho.T = 0.18, rho.C = 0.18),
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
    },
    formula = treatment ~ tte(eventtime, status, threshold = 2) + cont(toxicity.num, threshold = 0.5),
    sample.size = c(90,91,92,93,94,95),
    n.rep = 10, seed = 10, cpus = 1,
    method.inference = "u-statistic", scoring.rule = "Gehan")
)

summary(e.power)
 ## > statistic   : net benefit (null hypothesis Delta=0)
 ##     endpoint threshold n.T n.C mean.estimate sd.estimate mean.se rejection.rate
 ## toxicity.num       0.5  90  90        0.2468       0.084  0.0852          0.801
 ##                         91  91        0.2473      0.0836  0.0847          0.798
 ##                         92  92        0.2474      0.0832  0.0843          0.809
 ##                         93  93        0.2472      0.0826  0.0838           0.81
 ##                         94  94        0.2475      0.0816  0.0834          0.827
 ##                         95  95        0.2474       0.081  0.0829           0.84

