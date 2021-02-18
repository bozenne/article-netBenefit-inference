# path <- "p:/Cluster/GPC/Article-inference-Ustatistic/"
# setwd(path)
options(width = 100)

library(riskRegression)
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
gg.scatter <- gg.scatter + geom_point()
gg.scatter <- gg.scatter + labs(colour="OC group")
gg.scatter <- gg.scatter + scale_shape_manual("Left-censored",
                                              values = c(0,19),
                                              labels = c("yes","no"))
gg.scatter <- gg.scatter + geom_point() + xlab("Patient id") + ylab("Estradiol (nmol/L)")
gg.scatter

ggsave(gg.scatter, file = file.path("Results","ILLUSTRATION-scatterplot.pdf"))

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
head(dt.estra)

# > results
# endpoint threshold total(%) favorable(%) unfavorable(%) neutral(%) uninf(%)   delta   Delta
# estra      0.05      100         6.34          72.73       7.44     13.5 -0.6639 -0.6639
# p.value   
# 0.001 **
   
## e.BT_GehanR <- BuyseTest(group ~ tte(estra2, status = "observed", threshold = 0.05, censoring = "right"),
##                          data = dt.estra, scoring.rule = "Gehan")
## summary(e.BT_GehanR)

## * answer to reviewer 1
e.cox <- coxph(Surv(max(estra)-estra, observed) ~ strata(group),
               data = dt.estra, x = TRUE, y = TRUE)
e.pred <- predictCox(e.cox, keep.newdata = TRUE)
gg <- autoplot(e.pred, type = "survival", plot = FALSE)
gg2 <- gg$plot + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25),
                                    labels = -c(0,0.25,0.5,0.75,1,1.25)+max(dt.estra$estra))
gg2 <- gg2 + ylab("Survival") + xlab("Estradiol (nmol/L)") + labs(colour="OC group", shape = "Left censored")
gg2 <- gg2 + scale_shape_manual(breaks = c(0,1), values = c(3,18), labels = c("yes","no"))
gg2 <- gg2 + theme(text = element_text(size=20))
ggsave(gg2, filename = "figures/survivalCurve-biomarker.pdf")




