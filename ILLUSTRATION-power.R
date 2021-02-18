### ILLUSTRATION-power.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 22 2020 (10:17) 
## Version: 
## Last-Updated: okt  2 2020 (18:31) 
##           By: Brice Ozenne
##     Update #: 51
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
dt.prodige <- as.data.table(read.csv(file.path(path.data,"prodige.csv"), sep = ";", header = TRUE))

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
dt.prodige[, strata := interaction(oms_r,loca_r)]
  
rm.cols <- setdiff(names(dt.prodige),c("OS","bras","toxicity","etat","strata"))
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
## * re-analysis
## ** point estimate
system.time(
    e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 0.5, operator = "<0"),
                      data = dt.prodige, method.inference = "none")
)
   ## user  system elapsed 
   ## 0.03    0.00    0.03 

system.time(
    e.BT_Ustat <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 0.5, operator = "<0") + strata,
                            data = dt.prodige[strata %in% c("0.1","1.1","0.2","1.2")], method.inference = "u-statistic")
)
   ## user  system elapsed 
   ## 0.10    0.00    0.11 

summary(e.BT_Ustat)
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]    p.value    
##       OS       2.0   100.00        56.88          26.49      16.61     0.02 0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## toxicity       0.5    16.63         6.53           4.41       5.69     0.00 0.0212 0.3251   [0.2017;0.4383] 6.4399e-07 ***
system.time(
    e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
                           data = dt.prodige, method.inference = "bootstrap", n.resampling = 1e4,
                           trace = 1)
)
##   user  system elapsed 
## 220.60    0.47  225.61 

summary(e.BT_boot)
 ##    endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  delta  Delta CI [2.5% ; 97.5%]    p.value    
 ##           OS       2.0   100.00        56.88          26.49      16.61     0.02 0.3039 0.3039   [0.1816;0.4204] < 2.22e-16 ***
 ## toxicity.num       0.5    16.63         6.59           6.08       3.96     0.00 0.0051 0.3090   [0.1865;0.4273] < 2.22e-16 ***

library(microbenchmark)
e.bench <- microbenchmark(none = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
                                           data = dt.prodige, method.inference = "none"),
                          Ustat = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
                                            data = dt.prodige, method.inference = "u-statistic"),
                          bootstrap = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
                                                data = dt.prodige, method.inference = "bootstrap"),
                          times = 5)
summary(e.bench, unit = "s")
## >        expr         min          lq        mean      median          uq         max neval cld
## 1      none  0.02995687  0.03716029  0.03944084  0.04019474  0.04189177  0.04800055     5  a 
## 2     Ustat  0.09199358  0.09573385  0.10499052  0.09789393  0.11656605  0.12276522     5  a 
## 3 bootstrap 23.13373495 23.19421713 24.29345613 23.53133771 25.09936082 26.50863005     5   b

if(FALSE){
    ## install.packages("http://cran.r-project.org/src/contrib/Archive/BuyseTest/BuyseTest_1.0.tar.gz", lib = "C:/Users/hpl802/Downloads/LIBRTEMPO", type = "source", repos = NULL)
    library(BuyseTest, lib.loc="C:/Users/hpl802/Downloads/LIBRTEMPO") ## v2
    e.BT <- BuyseTest(treatment = "bras", endpoint = c("OS","toxicity"), type = c("TTE","cont"), censoring = c("etat",NA), threshold = c(2,0.5),
                      method = "Peron", data = dt.prodige[strata%in%c("0.1","1.1","0.2","1.2")], n.bootstrap = 0, strata = "strata")
    e.BT
    dt.prodige[, strata := interaction(oms_r,loca_r)]
    dt.prodige[,table(strata,bras)]
}

## ** sensitivity analysis
seqTh <- seq(0,6,0.5)
n.seqTh <- length(seqTh)
n.data <- NROW(dt.prodige)

ls.e.BT <- lapply(seqTh, function(iTh){
    ff <- eval(parse(text = paste0("bras ~ tte(OS, status = etat, threshold = ",iTh,")  + cont(toxicity.num, threshold = 0.5)")))
    BuyseTest(as.formula(ff),
              trace = 0,
              data = dt.prodige, method.inference = "u-statistic")
})
A.iid <- array(NA, dim = c(n.data, n.seqTh, 2))
A.iid[,,1] <- do.call(cbind,lapply(ls.e.BT,function(iObj){
    getIid(iObj, endpoint = 1)[[1]][,"favorable"]-getIid(iObj, endpoint = 1)[[1]][,"unfavorable"]
}))
A.iid[,,2] <- do.call(cbind,lapply(ls.e.BT,function(iObj){
    getIid(iObj, endpoint = 2)[[1]][,"favorable"]-getIid(iObj, endpoint = 2)[[1]][,"unfavorable"]
}))
M.beta <- do.call(cbind,lapply(ls.e.BT,function(iObj){
    coef(iObj)
}))
M.se <- do.call(cbind,lapply(ls.e.BT,function(iObj){
    sqrt(iObj@covariance[,"netBenefit"])
}))
sqrt(diag(crossprod(A.iid[,,1])))
sqrt(diag(crossprod(A.iid[,,2])))

dt.res <- as.data.table(transformCIBP(estimate = M.beta,
                        se = M.se,
                        iid = A.iid,
                        null = 0,
                        conf.level = 0.95,
                        alternative = "two.sided",
                        ci = TRUE, type = "none", min.value = -1, max.value = 1,
                        band = 1, method.band = "maxT-simulation", n.sim = 1e3, seed = 10,
                        p.value = FALSE),
                        estimate = M.beta,
                        se = M.se)
dt.res[["estimate"]] <- M.beta
dt.res[["se"]] <- M.se
dt.res[,c("row") := c("survival","survival+toxicity")[.SD$row]]
dt.res[,c("time") := seqTh[.SD$time]]
setnames(dt.res, old = c("row","time"), new = c("endpoint","threshold"))

gg <- ggplot(dt.res, aes(x=threshold))
gg <- gg + geom_line(aes(y=estimate, group = endpoint, color = endpoint, linetype = "estimate")) + geom_point(aes(y=estimate, shape = endpoint, color = endpoint))
## gg <- gg + geom_line(aes(y=lower, group = endpoint, linetype = "ci")) + geom_line(aes(y=upper, group = endpoint, linetype = "ci"))
gg <- gg + geom_line(aes(y=lowerBand, group = endpoint, color = endpoint, linetype = "band")) + geom_line(aes(y=upperBand, group = endpoint, color = endpoint, linetype = "band"))
gg <- gg + geom_hline(aes(yintercept = 0))
gg <- gg + xlab("Threshold of minimal clinical significance for survival (months)")
gg <- gg + ylab("Net benefit")
gg <- gg + scale_linetype_manual(name = NULL,
                                labels = c("simultaneous confidence intervals","point estimate"),
                                values = c(2,1))
gg <- gg + theme(text = element_text(size=10))

coef(ls.e.BT[[1]], statistic = "netBenefit")
confint(ls.e.BT[[1]])


lapply(ls.e.BT,coef)

## * simulation parameters
AFT0 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == 0,], dist = "Weibull")
AFT1 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == 1,], dist = "Weibull")
param0 <- exp(coef(AFT0))
param1 <- exp(coef(AFT1))

AFT0.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == 0,], dist = "Weibull")
AFT1.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == 1,], dist = "Weibull")
param0.cens <- exp(coef(AFT0.cens))
param1.cens <- exp(coef(AFT1.cens))


e.coxph <- coxph(Surv(OS,etat)~strata(bras), data = dt.prodige, x = TRUE, y = TRUE)
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
ggsave(gg.coxph, filename = file.path("p:/Cluster/GPC/Article-inference-Ustatistic/Results/ILLUSTRATION-folfirinox-surv.pdt"), height = 5, width = 7.5)


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
