## * data + package
source("ILLUSTRATION-data-management.R")

## * check values from the article
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
## censoring
dt.prodige[,100*mean(etat==0), by = "bras"]

## ustat
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 0.5, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic")
summary(e.BT)
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1583;0.3982] 1.3650e-05 ***

## boostrap
e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 0.5, operator = "<0"),
                       data = dt.prodige, method.inference = "bootstrap", n.resampling = 1e4, seed = 10,
                       trace = 1)
summary(e.BT_boot)
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
 ##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1852;0.4206] < 2.22e-16 ***
 ## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1607;0.4018] < 2.22e-16 ***

## * timing
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

## * sensitivity analysis
## several thresholds for the survival
seqTh <- seq(0,6,0.5)
n.seqTh <- length(seqTh)
n.data <- NROW(dt.prodige)

ls.e.BT <- lapply(seqTh, function(iTh){
    ff <- eval(parse(text = paste0("bras ~ tte(OS, status = etat, threshold = ",iTh,")  + cont(toxicity, threshold = 0.5, operator = \"<0\")")))
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

dt.res[,c("endpoint") := c("survival","survival+toxicity")[.SD$row]]
dt.res[,c("threshold") := seqTh[.SD$time]]

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

## * [not used] stratified analysis
if(FALSE){
    e.BT_strata <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 0.5, operator = "<0") + strata,
                             data = dt.prodige[strata %in% c("0.1","1.1","0.2","1.2")],
                             method.inference = "u-statistic")
    summary(e.BT_strata)
    ## endpoint threshold strata total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]   p.value   
    ##       OS         2 global   100.00        54.10          27.15      18.60     0.15  0.2695 0.2695   [0.1101;0.4153] 0.0010880 **
    ##                       0.1     7.59         3.99           2.20       1.33     0.07  0.2351                                      
    ##                       1.1    14.75         9.39           3.50       1.87     0.00  0.3993                                      
    ##                       0.2    17.07         9.77           4.22       3.01     0.07  0.3248                                      
    ##                       1.2    60.58        30.96          17.23      12.39     0.00  0.2266                                      
    ## toxicity       0.5 global    18.75         4.81           7.47       6.47     0.00 -0.0266 0.2429   [0.0823;0.3911] 0.0033069 **
    ##                       0.1     1.40         0.49           0.53       0.39     0.00 -0.0055                                      
    ##                       1.1     1.87         0.64           0.57       0.66     0.00  0.0049                                      
    ##                       0.2     3.08         0.48           1.58       1.03     0.00 -0.0643                                      
    ##                       1.2    12.39         3.20           4.79       4.40     0.00 -0.0263                                      
}
