# path <- "p:/Cluster/GPC/Article-inference-Ustatistic/"
setwd(path)
source("FCT-Bebu-Lachin.R")

library(BuyseTest)
library(pbapply)
library(data.table)
library(xtable)

## * example
set.seed(10)
d <- simBuyseTest(5)

## ** U-statistic
e.BT <- BuyseTest(treatment ~ cont(score), data = d, method.inference = "u-statistic", trace = FALSE)
out.BT <- list(estimate = c(favorable = as.double(coef(e.BT, statistic = c("favorable"))),
                            unfavorable = as.double(coef(e.BT, statistic = c("unfavorable")))),
               variance = matrix(c(e.BT@covariance[1,"favorable"],e.BT@covariance[1,"covariance"],e.BT@covariance[1,"covariance"],e.BT@covariance[1,"unfavorable"]),
                                 nrow = 2, ncol = 2,
                                 dimnames = list(c("favorable","unfavorable"),
                                                 c("favorable","unfavorable"))
                                 ))

## ** Manual
out.original <- UstatBL(formula = treatment ~ cont(score), data = d, method = "original") ## original formula of the article Bebu et al. (2015)
out.alternative <- UstatBL(formula = treatment ~ score, data = d, method = "alternative") ## another implementation of the formula of the article Bebu et al. (2015)
out.corrected <- UstatBL(formula = treatment ~ score, data = d, method = "corrected") ##  fix for the other implementation of the article Bebu et al. (2015)


## * simulation

warper <- function(i,n){ ## n <- 10
    set.seed(i)
    d <- simBuyseTest(n)
    e.BT <- BuyseTest(treatment ~ cont(score), data = d, method.inference = "u-statistic", trace = FALSE)
    out.original <- UstatBL(formula = treatment ~ cont(score), data = d, method = "original")
    out <- c(i = i, n = n,
             BT1 = confint(e.BT, order.Hprojection = 1)[1,c("estimate","se")],
             BT2 = confint(e.BT, order.Hprojection = 2)[1,c("estimate","se")],
             bebu = c(estimate = as.double(out.original$estimate["netBenefit"]), se = as.double(out.original$sd["netBenefit"])))
    return(out)
}

n.sim <- 5000
BuyseTest.options(order.Hprojection = 2)
ls.out <- c(pblapply(1:n.sim, warper, n = 5),
            pblapply(1:n.sim, warper, n = 10),
            pblapply(1:n.sim, warper, n = 20),
            pblapply(1:n.sim, warper, n = 30))
## |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 13m 01s
## |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 12m 13s
## |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 12m 17s
## |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 13m 20s
## saveRDS(ls.out, file = file.path("Results","SM-comparison.rds"))
dt.out <- as.data.table(do.call(rbind, ls.out))


## mean(is.na(dt.out$bebu.se))


tableSM1 <- dt.out[, .(empirical = sd(BT1.estimate),
                     BT1.se = mean(BT1.se),
                     BT2.se = mean(BT2.se),
                     bebu.se = mean(bebu.se, na.rm=TRUE)),
                 by = "n"]
names(tableSM1) <- c("sample size", "empirical","H-projection (order 1)","H-projection (order 2)","Bebu and Lachin")
tableSM1
##    sample size empirical H-projection (order 1) H-projection (order 2) Bebu and Lachin
## 1:           5 0.3781066              0.3534719              0.3656386       0.2684998
## 2:          10 0.2678378              0.2561108              0.2613307       0.2286741
## 3:          20 0.1841571              0.1822960              0.1843487       0.1729108
## 4:          30 0.1503660              0.1489467              0.1501046       0.1438876
xtable(tableSM1, digits = 4)

dt.out[,BT1.pvalue := 2*(1-pnorm(abs(BT1.estimate)/BT1.se))]
dt.out[,BT2.pvalue := 2*(1-pnorm(abs(BT2.estimate)/BT2.se))]
dt.out[,bebu.pvalue := 2*(1-pnorm(abs(bebu.estimate)/bebu.se))]
dt.out[,.(BT1.type1 = mean(BT1.pvalue<=0.05),
          BT2.type1 = mean(BT2.pvalue<=0.05),
          bebu.type1 = mean(bebu.pvalue<=0.05, na.rm = TRUE)), by = "n"]
##     n BT1.type1 BT2.type1 bebu.type1
## 1:  5    0.1326    0.1196  0.1836653
## 2: 10    0.0880    0.0824  0.1208000
## 3: 20    0.0660    0.0630  0.0804000
## 4: 30    0.0620    0.0598  0.0698000

## * Example where the formula of Bebu et al. (2015) gives a negative variance
d <- data.table(treatment = c(rep("C",5),rep("T",5)),
                score = c(9,7,8,10,5,4,6,3,1,2))
out.original <- UstatBL(formula = treatment ~ cont(score), data = d, method = "original") ## original formula of the article Bebu et al. (2015)
BuyseTest.options(order.Hprojection = 2)
e.BT <- BuyseTest(treatment ~ cont(score),
                  data = d,
                  method.inference = "u-statistic",
                  trace = FALSE)
confint(e.BT, statistic = "favorable", order.Hprojection = 2)[,"se"]^2
confint(e.BT, statistic = "favorable", order.Hprojection = 1)[,"se"]^2

## manual calculation
X <- c(4,6,3,1,2)
Y <- c(9,7,8,10,5)
m <- 5
n <- 5
N <- m+n

grid <- cbind(expand.grid(i = 1:m, j = 1:n),expand.grid(Xi = X,Yj = Y))
grid$Sij1 <- apply(grid,1, function(x){x["Xi"]>x["Yj"]})
grid$Sij2 <- apply(grid,1, function(x){x["Xi"]<x["Yj"]})
dt.grid <- as.data.table(grid)

Si <- dt.grid[,.(Si1 = sum(Sij1), Si2 = sum(Sij2)),by = "i"]
pXY1 <- sum(dt.grid$Sij1) / (m*n) ## P[X1>Y1] = 1/mn \sum_ij S1_ij

## P[X1>Y1,X1>Y1']
sum(Si$Si1*(Si$Si1-1)) / (m*n*(n-1)) ## = 0
gridJ <- cbind(expand.grid(i = 1:m, j1 = 1:n, j2 = 1:n), expand.grid(Xi = X, Yj1 = Y, Yj2 = Y))
mean(apply(gridJ,1, function(x){(x["Xi"]>x["Yj1"])*(x["Xi"]>x["Yj2"])})) ## > 0



######################################################################
### SM-comparison-bebu-lachin.R ends here
