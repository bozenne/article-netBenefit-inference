## * data + package
library(exact2x2)
source("ILLUSTRATION-0-data-management.R")
options(width = 120)

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
round(100*prop.table(dt.prodige[,table(toxicity, bras)],2),2)
##         bras
## toxicity Gemcitabine Folfirinox
##        0        1.17       3.51
##        1        2.92       4.09
##        2       36.26      23.39
##        3       39.18      47.37
##        4       19.88      21.05
##        5        0.58       0.58

colSums(100*prop.table(dt.prodige[,table(toxicity, bras)],2)[4:6,])
## Gemcitabine  Folfirinox 
##    59.64912    69.00585 

uncondExact2x2(x1 = table(dt.prodige$toxicity>=3,dt.prodige$bras)["TRUE","Gemcitabine"],
               n1 = sum(dt.prodige$bras=="Gemcitabine"),
               x2 = table(dt.prodige$toxicity>=3,dt.prodige$bras)["TRUE","Folfirinox"],
               n2 = sum(dt.prodige$bras=="Folfirinox"), conf.int = TRUE)
## data:  x1/n1=(102/171) and x2/n2= (118/171)
## proportion 1 = 0.59649, proportion 2 = 0.69006, p-value = 0.07418
## alternative hypothesis: true p2-p1 is not equal to 0
## 95 percent confidence interval:
##  -0.009435654  0.200794220
## sample estimates:
##      p2-p1 
## 0.09356725 


fisher.test(dt.prodige[,table(toxicity, bras)])

## * re-analysis
## censoring
dt.prodige[,.(censoring=100*mean(etat==0)), by = "bras"]
##           bras censoring
## 1: Gemcitabine  14.03509
## 2:  Folfirinox  26.31579

## ustat
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 1, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic")
summary(e.BT)
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1819;0.4166] 2.1557e-06 ***
## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1583;0.3982] 1.3650e-05 ***



## boostrap
e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 0.5, operator = "<0"),
                       data = dt.prodige, method.inference = "bootstrap", n.resampling = 1e4, seed = 10, cpus = 10,
                       trace = 1)
summary(e.BT_boot)
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta  Delta CI [2.5% ; 97.5%]    p.value    
 ##       OS       2.0   100.00        56.88          26.49      16.61     0.02  0.3039 0.3039   [0.1852;0.4206] < 2.22e-16 ***
 ## toxicity       0.5    16.63         4.41           6.53       5.69     0.00 -0.0212 0.2827   [0.1607;0.4018] < 2.22e-16 ***

## 6.5/16.63
## [1] 0.3908599
## > 4.41/16.61
## [1] 0.2655027

confint(e.BT)
##             estimate         se  lower.ci  upper.ci      p.value
## OS_2       0.3038845 0.06014829 0.1818654 0.4166871 2.188018e-06
## toxicity_1 0.2826621 0.06152741 0.1581706 0.3983080 1.390671e-05
confint(e.BT_boot)
##               estimate         se  lower.ci  upper.ci p.value
## OS_2         0.3038845 0.06076773 0.1825247 0.4221415       0
## toxicity_0.5 0.2826621 0.06218618 0.1571735 0.4034210       0

## * timing
## library(microbenchmark)
## e.bench <- microbenchmark(none = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
##                                            data = dt.prodige, method.inference = "none"),
##                           Ustat = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
##                                             data = dt.prodige, method.inference = "u-statistic"),
##                           bootstrap = BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity.num, threshold = 0.5),
##                                                 data = dt.prodige, method.inference = "bootstrap"),
##                           times = 5)
## summary(e.bench, unit = "s")
## ## >        expr         min          lq        mean      median          uq         max neval cld
## ## 1      none  0.02995687  0.03716029  0.03944084  0.04019474  0.04189177  0.04800055     5  a 
## ## 2     Ustat  0.09199358  0.09573385  0.10499052  0.09789393  0.11656605  0.12276522     5  a 
## ## 3 bootstrap 23.13373495 23.19421713 24.29345613 23.53133771 25.09936082 26.50863005     5   b

system.time(
    BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 0.5, operator = "<0"),
              data = dt.prodige, method.inference = "bootstrap", n.resampling = 1e4, seed = 10, cpus = 1,
              trace = 1)
)
##    user  system elapsed 
## 148.399   0.039 148.479 

microbenchmark::microbenchmark(
                    BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 1, operator = "<0"),
                              data = dt.prodige, method.inference = "u-statistic", trace = FALSE)
                )
## Unit: milliseconds
##                                                                                                                                                                                    expr
##  BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity,      threshold = 1, operator = "<0"), data = dt.prodige, method.inference = "u-statistic",      trace = FALSE)
##       min       lq     mean   median       uq      max neval
##  41.12059 41.82398 44.99263 44.25925 45.45983 133.0568   100

## * sensitivity analysis
system.time(
    eSe.BT <- sensitivity(e.BT, threshold = list(OS=seq(0,12,1), toxicity = 1:3), band = FALSE, adj.p.value = FALSE)
)
##  user  system elapsed 
## 1.866   0.000   1.866 
system.time(
    eSe.BT <- sensitivity(e.BT, threshold = list(OS=seq(0,12,1), toxicity = 1:3), band = TRUE, adj.p.value = FALSE)
)
## user  system elapsed 
## 3.203   0.012   3.229
## autoplot(eSe.BT)

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
