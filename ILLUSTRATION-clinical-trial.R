### ILLUSTRATION-clinical-trial.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 22 2020 (10:17) 
## Version: 
## Last-Updated: okt  8 2020 (17:50) 
##           By: Brice Ozenne
##     Update #: 53
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

## * Data management
## ** load
dt.prodige <- as.data.table(read.csv(file.path(path.data,"PRODIGE.csv"), sep = ";", header = TRUE))
## dt.prodige <- as.data.table(read.csv("data/PRODIGE.csv", sep = ";", header = TRUE))

## ** process
dt.prodige[, d_dn2 := as.Date(d_dn, "%d/%m/%Y")]
dt.prodige[, randodt2 := as.Date(randodt, "%d/%m/%Y")]
dt.prodige[, bras := factor(ifelse(bras==2,"Gemcitabine","Folfirinox"), c("Gemcitabine","Folfirinox"))]
dt.prodige[, OS := as.numeric(difftime(d_dn2,randodt2,units="days")/30.44)]

M.tox <- cbind(inf = dt.prodige[,pmax(mx_inf1, mx_inf2, mx_inf3, mx_inf4, 0, na.rm = TRUE)],
               car = dt.prodige[,pmax(mx_car1, mx_car2, mx_car3, mx_car4, mx_car5, mx_car6, mx_car7, 0, na.rm = TRUE)],
               aon = dt.prodige[,pmax(mx_aon1, mx_aon2, mx_aon3, mx_aon4, mx_aon5, mx_aon6, mx_aon7, 0, na.rm = TRUE)],
               pul = dt.prodige[,pmax(mx_pul1, mx_pul2, mx_pul3, mx_pul4, 0, na.rm = TRUE)],
               aut = dt.prodige[,pmax(mx_aut1, mx_aut2, mx_aut3, mx_aut4, 0, na.rm = TRUE)],
               ren = dt.prodige[,pmax(mx_ren1, mx_ren2, mx_ren3, mx_ren4, 0, na.rm = TRUE)],
               gas = dt.prodige[,pmax(mx_gas1, mx_gas2, mx_gas3, mx_gas4, mx_gas5, mx_gas6, mx_gas7, 0, na.rm = TRUE)],
               der = dt.prodige[,pmax(mx_der1, mx_der2, mx_der3, mx_der4, mx_der5, mx_der6, mx_der7, 0, na.rm = TRUE)],
               tog = dt.prodige[,pmax(mx_tog1, mx_tog2, mx_tog3, mx_tog4, 0, na.rm = TRUE)],
               uro = dt.prodige[,pmax(mx_uro1, mx_uro2, mx_uro3, mx_uro4, 0, na.rm = TRUE)]
               )
dt.prodige[,toxicity := apply(M.tox,1, max)]
dt.prodige[, strata := interaction(oms_r,loca_r)] ## severity of illness (Eastern Cooperative Oncology Group) and primary tumor localization (the head vs. the body or tail of the pancreas)
  
rm.cols <- setdiff(names(dt.prodige),c("OS","bras","toxicity","etat","strata","oms_r","loca_r"))
dt.prodige[,c(rm.cols) := NULL]


## ** check values
dt.prodige[,.N, by = "bras"]
##           bras   N
## 1: Gemcitabine 171
## 2:  Folfirinox 171


e.coxph <- coxph(Surv(OS, etat)~bras, data = dt.prodige)
summary(e.coxph)
##                exp(coef) exp(-coef) lower .95 upper .95
## brasFolfirinox    0.5674      1.762     0.446    0.7219

e.coxphS <- coxph(Surv(OS, etat)~strata(bras), data = dt.prodige, x = TRUE, y = TRUE)
predictCox(e.coxphS, times = c(11.1,6.7))
   ## observation      strata times cumhazard survival
## 1:           1 Gemcitabine   6.7     0.669    0.512
## 4:           4  Folfirinox  11.1     0.676    0.509

dt.prodige[,table(toxicity, bras)]
##         bras
## toxicity Gemcitabine Folfirinox
##        0           2          6
##        1           5          7
##        2          62         40
##        3          67         81
##        4          34         36
##        5           1          1


## * Parametric model
AFT_G <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[bras == "Gemcitabine",], dist = "Weibull")
AFT_F <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[bras == "Folfirinox",], dist = "Weibull")
param_G <- exp(coef(AFT_G))
param_F <- exp(coef(AFT_F))

AFT.cens_G <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[bras == "Gemcitabine",], dist = "Weibull")
AFT.cens_F <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[bras == "Folfirinox",], dist = "Weibull")
param.cens_G <- exp(coef(AFT.cens_F))
param.cens_F <- exp(coef(AFT.cens_G))


## * Figure
e.coxph <- coxph(Surv(OS,etat)~strata(bras), data = dt.prodige, x = TRUE, y = TRUE)
pred.coxph <- predictCox(e.coxph, type = "survival", keep.newdata = TRUE)
gg.coxph <- autoplot(pred.coxph,plot = FALSE)$plot


X <- seq(0,50, by = 0.1)
dt <- rbind(data.table(strata = "Gemcitabine",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param_G["scale"], shape = param_G["shape"])
                       ),
            data.table(strata = "Folfirinox",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param_F["scale"], shape = param_F["shape"])
                       )
            )

gg.coxph <- gg.coxph + geom_line(data = dt, mapping = aes(x = time, y = surv, group = strata, color = type), size = 0.8)
gg.coxph$layers <- gg.coxph$layers[c(3,1,2)]
gg.coxph <- gg.coxph + scale_colour_manual(name = "",
                                           values = c("Weibull model" = "black",
                                                      "Gemcitabine" = "darkorange",
                                                      "Folfirinox" = "darkblue"
                                                      ),
                                           labels = c("Weibull model" = "Weibull model",
                                                      "Gemcitabine" = "arm Gemcitabine",
                                                      "Folfirinox" = "arm Folfirinox"
                                                      ))
gg.coxph <- gg.coxph + scale_shape_manual(name = "",
                                          breaks = c(0,1),
                                          values = c(3,18),
                                          labels = c("censoring","death"))
gg.coxph <- gg.coxph + theme(legend.position="right", text = element_text(size=15))
gg.coxph <- gg.coxph + xlab("Months") + ylab("Probability of survival")
ggsave(gg.coxph, filename = file.path("p:/Cluster/GPC/Article-inference-Ustatistic/Results/ILLUSTRATION-folfirinox-surv.pdf"), height = 5, width = 7.5)


## * Power calculation
## **  estimation of rho
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

    dt.sim <- simBuyseTest(n,
                           argsBin = argsBin,
                           argsTTE = argsTTE,
                           latent = TRUE)
    dt.sim$toxicity.num <- as.numeric(dt.sim$toxicity)-1
    dt.sim$one <- 1

    ## Check marginal distribution
    ## plot(prodlim(Hist(OS, etat) ~ bras, data = dt.prodige), legend = FALSE, atrisk = FALSE)
    ## plot(prodlim(Hist(eventtime, status) ~ treatment, data = dt.sim), legend = FALSE, atrisk = FALSE, add = TRUE)

    ## dt.sim[treatment == "C", table(toxicity)/.N]
    ## dt.sim[treatment == "T", table(toxicity)/.N]
    iE.BT <- BuyseTest(treatment ~ tte(eventtimeUncensored, status = one, threshold = 2) + cont(toxicity.num, threshold = 0.5),
                       data = dt.sim)
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

## ** simulation
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


## *** 1 CPUS, n.rep = 100
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

## *** 4 CPUS, n.rep = 1000
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

## * Benefit Risk analysis 
system.time(
    e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, operator = "<0") + strata,
                      data = dt.prodige[oms_r %in% 0:1], method.inference = "none")
)
   ## user  system elapsed 
   ## 0.05    0.00    0.05 
summary(e.BT, strata = "global")
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta
 ##       OS     2e+00   100.00        54.10          27.15      18.60     0.15  0.2695 0.2695
 ## toxicity     1e-12    18.75         4.81           7.47       6.47     0.00 -0.0266 0.2429

system.time(
    e.BT_Ustat <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, operator = "<0") + strata,
                            data = dt.prodige[oms_r %in% 0:1], method.inference = "u-statistic")
)
   ## user  system elapsed 
   ##  0.1     0.0     0.1 
summary(e.BT_Ustat, strata = "global")
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]   p.value   
##       OS     2e+00   100.00        54.10          27.15      18.60     0.15  0.2695 0.2695   [0.1101;0.4153] 0.0010880 **
## toxicity     1e-12    18.75         4.81           7.47       6.47     0.00 -0.0266 0.2429   [0.0823;0.3911] 0.0032956 **
system.time(
    e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, operator = "<0") + strata,
                           data = dt.prodige[oms_r %in% 0:1], method.inference = "bootstrap", strata.resampling = "strata", n.resampling = 1e4,
                           trace = 1, cpus = 3, seed = 10)
)
   ## user  system elapsed 
   ## 6.05    1.86  221.38 
summary(e.BT_boot, strata = "global")
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]   p.value    
 ##       OS     2e+00   100.00        54.10          27.15      18.60     0.15  0.2695 0.2695   [0.1373;0.3954]  <2e-16 ***
 ## toxicity     1e-12    18.75         4.81           7.47       6.47     0.00 -0.0266 0.2429   [0.1073;0.3728]   0.001  **

system.time(
    e.BT_perm <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, operator = "<0") + strata,
                           data = dt.prodige[oms_r %in% 0:1], method.inference = "studentized permutation", strata.resampling = "strata", n.resampling = 1e4,
                           trace = 1, cpus = 3, seed = 10)
)
   ## user  system elapsed 
   ## 5.96    3.04  366.79 

summary(e.BT_perm, strata = "global")
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta    p.value    
 ##       OS     2e+00   100.00        54.10          27.15      18.60     0.15  0.2695 0.2695 < 2.22e-16 ***
 ## toxicity     1e-12    18.75         4.81           7.47       6.47     0.00 -0.0266 0.2429 < 2.22e-16 ***

## rbind(confint(e.BT_Ustat),
##       confint(e.BT_boot),
##       confint(e.BT_perm))
      


######################################################################
### ILLUSTRATION-clinical-trial.R ends here
