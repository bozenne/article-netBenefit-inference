# path <- "p:/Cluster/GPC/Article-inference-Ustatistic/"
# setwd(path)
options(width = 100)

library(riskRegression)
library(prodlim)
library(survival)
library(data.table)
library(BuyseTest)
library(ggplot2)

## * Data management
df.all <- read.csv(file.path("data","Final.OC.csv"), sep=";", dec=",", na.strings = "")
df.estra <- df.all[!is.na(df.all$Imputed.ESTRA),c("CIMBI.ID","OC.1yes","LL.ESTRA","Imputed.ESTRA")]
dt.estra <- as.data.table(df.estra)
dt.estra[, OC.1yes := factor(OC.1yes, levels = 0:1, labels = c("non-user","user"))]
setnames(dt.estra, old = "LL.ESTRA", new = "censored")
setnames(dt.estra, old = "Imputed.ESTRA", new = "estra")
setnames(dt.estra, old = "OC.1yes", new = "group")
setnames(dt.estra, old = "CIMBI.ID", new = "id")
dt.estra[, observed := 1-censored]
dt.estra[, observed.char := factor(observed, 0:1)]
dt.estra[, censored := NULL]
dt.estra[, estra2 := ceiling(max(estra))+1-estra]
## dt.estra[, estra2 := -estra]
setkeyv(dt.estra, cols = c("group","estra2"))
dt.estra[,id2 := 1:.N]

## * Descriptive statistics
dt.estra[, .(n = .N, observed = sum(observed), observedPC = 100*mean(observed)),by = "group"]
#      group  n observed observedPC
#1: non-user 33       28   84.84848
#2:     user 11        5   45.45455

gg.scatter <- ggplot(dt.estra,aes(x = id2, y = estra, color = group, shape = observed.char))
gg.scatter <- gg.scatter + geom_point(size = 3)
gg.scatter <- gg.scatter + labs(colour="OC group")
gg.scatter <- gg.scatter + scale_shape_manual("Censored",
                                              values = c(0,19),
                                              labels = c("yes","no"))
gg.scatter <- gg.scatter + geom_point() + xlab("Patient id") + ylab("Estradiol (nmol/L)")
gg.scatter <- gg.scatter + theme(text = element_text(size=25),
                                 axis.line = element_line(size = 1.25),
                                 axis.ticks = element_line(size = 2),
                                 axis.ticks.length=unit(.25, "cm"))
gg.scatter

if(FALSE){
    ggsave(gg.scatter, file = file.path("figures","ILLUSTRATION2-scatterplot.pdf"), width = 8, height = 5)
}

gg.hist <- ggplot(dt.estra,aes(estra, fill = group)) + geom_histogram(breaks = seq(0,1.25,by = 0.025))

dt.estra[, sd(estra),by = "group"]
# group         V1
# 1: non-user 0.31279108
# 2:     user 0.06121497
dt.estra[, mean(estra),by = "group"]
# group         V1
# 1: non-user 0.37151515
# 2:     user 0.09454545


## * GPC
BuyseTest.options(order.Hprojection = 2)
e.BT_GehanL <- BuyseTest(group ~ tte(estra, status = observed, threshold = 0.05, censoring = "left"),
                         data = dt.estra,
                         scoring.rule = "Gehan")
summary(e.BT_GehanL)
 ## - statistic       : net benefit (delta: endpoint specific, Delta: global) 
 ## - null hypothesis : Delta == 0 
 ## - confidence level: 0.95 
 ## - inference       : H-projection of order 2
 ## - treatment groups: non-user (control) vs. user (treatment) 
 ## - censored pairs  : deterministic score or uninformative
 ## - uninformative pairs: no contribution
 ## - results
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   Delta CI [2.5% ; 97.5%]    p.value    
 ##    estra      0.05      100         6.34          72.73       7.44     13.5 -0.6639 [-0.8255;-0.4018] 2.7728e-05 ***
confint(e.BT_GehanL)
##              estimate        se   lower.ci   upper.ci      p.value
## estra_0.05 -0.6639118 0.1067079 -0.8254762 -0.4017934 2.772777e-05

## same as
dt.estra[, estra2 := ceiling(max(estra))+1-estra]
e.BT_GehanR <- BuyseTest(group ~ tte(estra2, status = "observed", threshold = 0.04999999, operator = "<0"),
                         data = dt.estra, scoring.rule = "Gehan", keep.pairScore = TRUE)
summary(e.BT_GehanR)
 ## endpoint  threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   Delta CI [2.5% ; 97.5%]    p.value    
 ##   estra2 0.04999999      100         6.34          72.73       7.44     13.5 -0.6639 [-0.8255;-0.4018] 2.7728e-05 ***
e.BT_PeronR <- BuyseTest(group ~ tte(estra2, status = "observed", threshold = 0.04999999, operator = "<0"),
                         data = dt.estra, scoring.rule = "Peron", keep.pairScore = TRUE)
summary(e.BT_PeronR)
 ## endpoint  threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   Delta CI [2.5% ; 97.5%]    p.value    
 ##   estra2 0.04999999      100         8.54          74.27       7.55     9.64 -0.6573 [-0.8294;-0.3712] 0.00010517 ***

dt.estra[c(29,36)]



## * answer to reviewer 1
e.cox <- coxph(Surv(max(estra)-estra, observed) ~ strata(group),
               data = dt.estra, x = TRUE, y = TRUE)
e.pred <- predictCox(e.cox, keep.newdata = TRUE)
gg <- autoplot(e.pred, type = "survival", plot = FALSE)
gg2 <- gg$plot + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25),
                                    labels = -c(0,0.25,0.5,0.75,1,1.25)+max(dt.estra$estra))
gg2 <- gg2 + ylab("Survival") + xlab("Estradiol (nmol/L)") + labs(colour="OC group", shape = "Left censored")
gg2 <- gg2 + scale_shape_manual(breaks = c(0,1), values = c(3,18), labels = c("yes","no"))
gg2 <- gg2 + coord_cartesian(ylim = c(0,1))
gg2 <- gg2 + theme(text = element_text(size=20))
gg2 <- gg2 + theme(text = element_text(size=25),
                   axis.line = element_line(size = 1.25),
                   axis.ticks = element_line(size = 2),
                   axis.ticks.length=unit(.25, "cm"))
if(FALSE){
    ggsave(gg2, filename = "figures/ILLUSTRATION2-survCurve.pdf", width = 8, height = 5)
}



