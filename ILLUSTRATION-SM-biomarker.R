# path <- "p:/Cluster/GPC/Article-inference-Ustatistic/"
# setwd(path)
options(width = 100)

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

dt.estra[, sd(estra),by = "group"]
# group         V1
# 1: non-user 0.31279108
# 2:     user 0.06121497
dt.estra[, mean(estra),by = "group"]
# group         V1
# 1: non-user 0.37151515
# 2:     user 0.09454545

## * BuyseTest
BuyseTest.options(order.Hprojection = 2)
e.BT_GehanL <- BuyseTest(group ~ tte(estra, status = "observed", threshold = 0.05, censoring = "left"),
                         data = dt.estra,
                         scoring.rule = "Gehan")
summary(e.BT_GehanL)
## - results
## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   Delta
##    estra      0.05      100         6.34          72.73       7.44     13.5 -0.6639
## CI [2.5% ; 97.5%]    p.value    
## [-0.8255;-0.4018] 2.7728e-05 ***
confint(e.BT_GehanL)
##              estimate        se   lower.ci   upper.ci      p.value
## estra_0.05 -0.6639118 0.1067079 -0.8254762 -0.4017934 2.772777e-05

## same as:
dt.estra$time <- 1e-3+max(dt.estra$estra)-dt.estra$estra
e.BT_GehanR <- BuyseTest(group ~ tte(time, status = "observed", threshold = 0.05, censoring = "right"),
                         data = dt.estra,
                         scoring.rule = "Gehan")
summary(e.BT_GehanR)
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta
 ##     time      0.05      100        72.73           6.34       7.44     13.5 0.6639
 ## CI [2.5% ; 97.5%]    p.value    
 ##   [0.4018;0.8255] 2.7728e-05 ***

e.BT_PeronR <- BuyseTest(group ~ tte(time, status = "observed", threshold = 0.05, censoring = "right"),
                         data = dt.estra,
                         scoring.rule = "Peron")
summary(e.BT_PeronR)
 ## endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)  Delta
 ##     time      0.05      100        74.16           8.54       7.66     9.64 0.6562
 ## CI [2.5% ; 97.5%]    p.value    
 ##   [0.3706;0.8284] 0.00010414 ***




