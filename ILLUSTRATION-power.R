### ILLUSTRATION-power.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 22 2020 (10:17) 
## Version: 
## Last-Updated: jun  6 2020 (14:51) 
##           By: Brice Ozenne
##     Update #: 35
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

path <- "p:/Cluster/GPC/Article-inference-Ustatistic/"
setwd(path)
options(width = 130)

path.data <- "data"

## * Package
library(BuyseTest) ## butils.base:::sourcePackage("BuyseTest", c.code = TRUE, trace = TRUE)
library(survival)
library(prodlim)
library(flexsurv)
library(riskRegression)
library(ggplot2)
library(ggthemes)
library(data.table)

## * data management
## ** load
df.prodige <- read.csv(file.path(path.data,"data PRODIGE.csv"), sep = ";", header = TRUE)

## ** process
df.prodige$d_dn2 <- as.Date(df.prodige$d_dn, "%d/%m/%Y")
df.prodige$randodt2 <- as.Date(df.prodige$randodt, "%d/%m/%Y")
df.prodige$bras <- ifelse(df.prodige$bras==2,0,1)
df.prodige$OS <- as.numeric(difftime(df.prodige$d_dn2,df.prodige$randodt2,units="days")/30.44)

## * re-analysis
system.time(
    e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                      data = df.prodige, method.inference = "none")
)
## user  system elapsed 
## 0.14    0.00    0.22 

system.time(
    e.BT_Ustat <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                            data = df.prodige, method.inference = "u-statistic")
)

summary(e.BT_Ustat)
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]    p.value    
##       OS         2      100        56.88          26.49      16.61     0.02 0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
system.time(
    e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                           data = df.prodige, method.inference = "bootstrap", n.resampling = 1e3,
                           trace = 1)
)
summary(e.BT_boot)
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]    p.value    
##       OS         2      100        56.88          26.49      16.61     0.02 0.3039 0.3039   [0.1797;0.4101] < 2.22e-16 ***

library(microbenchmark)
e.bench <- microbenchmark(none = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                                           data = df.prodige, method.inference = "none"),
                          Ustat = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                                            data = df.prodige, method.inference = "u-statistic"),
                          bootstrap = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2),
                                                data = df.prodige, method.inference = "bootstrap"),
                          times = 5)
summary(e.bench, unit = "s")
##        expr         min         lq        mean     median          uq        max neval cld
## 1      none  0.03884476  0.0804832  0.07901163  0.0845405  0.08753235  0.1036573     5  a 
## 2     Ustat  0.28849105  0.2895154  0.35713805  0.2909330  0.32857495  0.5881758     5  a 
## 3 bootstrap 33.00206102 53.8771571 56.58843587 62.3667257 66.63153893 67.0646965     5   b

## * simulation parameters
AFT0 <- flexsurvreg(Surv(OS, etat) ~ 1, data = df.prodige[df.prodige$bras == 0,], dist = "Weibull")
AFT1 <- flexsurvreg(Surv(OS, etat) ~ 1, data = df.prodige[df.prodige$bras == 1,], dist = "Weibull")
param0 <- exp(coef(AFT0))
param1 <- exp(coef(AFT1))

AFT0.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = df.prodige[df.prodige$bras == 0,], dist = "Weibull")
AFT1.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = df.prodige[df.prodige$bras == 1,], dist = "Weibull")
param0.cens <- exp(coef(AFT0.cens))
param1.cens <- exp(coef(AFT1.cens))


e.coxph <- coxph(Surv(OS,etat)~strata(bras), data = df.prodige, x = TRUE, y = TRUE)
pred.coxph <- predictCox(e.coxph, type = "survival", keep.newdata = TRUE)
gg.coxph <- autoplot(pred.coxph,plot = FALSE)$plot

X <- seq(0,50, by = 0.1)
dt <- rbind(data.table(strata = "bras=0",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param0["scale"], shape = param0["shape"])
                       ),
            data.table(strata = "bras=1",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param1["scale"], shape = param1["shape"])
                       )
            )

gg.coxph <- gg.coxph + geom_line(data = dt, mapping = aes(x = time, y = surv, group = strata, color = type), size = 0.8)
gg.coxph$layers <- gg.coxph$layers[c(3,1,2)]
gg.coxph <- gg.coxph + scale_colour_manual(name = "",
                                           values = c("Weibull model" = "black",
                                                      "bras=0" = "darkorange",
                                                      "bras=1" = "darkblue"
                                                      ),
                                           labels = c("Weibull model" = "Weibull model",
                                                      "bras=0" = "arm Gemcitabine",
                                                      "bras=1" = "arm Folfirinox"
                                                      ))
gg.coxph <- gg.coxph + scale_shape_manual(name = "",
                                          breaks = c(0,1),
                                          values = c(3,18),
                                          labels = c("censoring","death"))
gg.coxph <- gg.coxph + theme(legend.position="right", text = element_text(size=15))
gg.coxph <- gg.coxph + xlab("Months") + ylab("Probability of survival")
ggsave(gg.coxph, filename = file.path("p:/Cluster/GPC/Article-inference-Ustatistic/Results/ILLUSTRATION-folfirinox-surv.pdf"), height = 5, width = 7.5)


## ** find rho
n <- 2e4
seqRho <- seq(0,0.3,by=0.03)

out <- NULL
for(iRho in 1:length(seqRho)){
    cat(iRho, " ")
    argsBin <- list(p.C = c("0" = 1.2/100,
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
                    rho.T = seqRho[iRho], rho.C = seqRho[iRho])

    argsTTE <- list(scale.C = 0.95*param0["scale"],
                    scale.T = 0.97*param1["scale"],
                    shape.C = param0["shape"],
                    shape.T = param1["shape"],
                    scale.Censoring.C = param0.cens["scale"],
                    scale.Censoring.T = param1.cens["scale"],
                    shape.Censoring.C = param0.cens["shape"],
                    shape.Censoring.T = param1.cens["shape"])

    df.sim <- simBuyseTest(n,
                           argsBin = argsBin,
                           argsTTE = argsTTE,
                           latent = TRUE)
    df.sim$toxicity.num <- as.numeric(df.sim$toxicity)-1
    df.sim$one <- 1

    ## Check marginal distribution
    ## plot(prodlim(Hist(OS, etat) ~ bras, data = df.prodige), legend = FALSE, atrisk = FALSE)
    ## plot(prodlim(Hist(eventtime, status) ~ treatment, data = df.sim), legend = FALSE, atrisk = FALSE, add = TRUE)

    ## df.sim[treatment == "C", table(toxicity)/.N]
    ## df.sim[treatment == "T", table(toxicity)/.N]
    iE.BT <- BuyseTest(treatment ~ tte(eventtimeUncensored, status = one, threshold = 2) + cont(toxicity.num, threshold = 0.5),
                       data = df.sim)
    out <- rbind(out,
                 cbind(rho = seqRho[iRho], delta1 = coef(iE.BT)[1], delta2 = coef(iE.BT)[2]-coef(iE.BT)[1], Delta = coef(iE.BT)[2])
                 )
}
rownames(out) <- NULL
out

##        rho    delta1      delta2     Delta
##  [1,] 0.00 0.2713407  0.01147245 0.2828132
##  [2,] 0.03 0.2661984 -0.00014414 0.2660543
##  [3,] 0.06 0.2681221 -0.00828872 0.2598334
##  [4,] 0.09 0.2767830 -0.01567913 0.2611039
##  [5,] 0.12 0.2649490 -0.02238958 0.2425594
##  [6,] 0.15 0.2633671 -0.02561992 0.2377471
##  [7,] 0.18 0.2737693 -0.02920194 0.2445673
##  [8,] 0.21 0.2816950 -0.03097143 0.2507236
##  [9,] 0.24 0.2831517 -0.03300977 0.2501420
## [10,] 0.27 0.2732482 -0.03662647 0.2366218
## [11,] 0.30 0.2713966 -0.03722738 0.2341692

## * simulation
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
)


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

######################################################################
### ILLUSTRATION-power.R ends here
