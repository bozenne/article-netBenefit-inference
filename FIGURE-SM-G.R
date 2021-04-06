## * path
## path <- "P:/Cluster/GPC/Article-inference-Ustatistic-Rao"
## setwd(path)
path.results <- "Results"
path.figures <- "figures-article"


## * libraries
library(survival)
library(riskRegression)
library(data.table)
library(ggplot2)
library(ggpubr) ## neededed for graphical display
library(ggthemes) ## neededed for graphical display
library(xtable) ## neededed for display
source("FCT-gg.R")


## * load data
source("ILLUSTRATION2-biomarker.R")

## * Create figure
## ** scatterplot
figureSMG.a <- ggplot(dt.estra,aes(x = id2, y = estra, color = group, shape = observed.char))
figureSMG.a <- figureSMG.a + geom_point(size = 3)
figureSMG.a <- figureSMG.a + labs(colour="OC group")
figureSMG.a <- figureSMG.a + scale_shape_manual("Below the detection limit",
                                              values = c(0,19),
                                              labels = c("yes","no"))
figureSMG.a <- figureSMG.a + geom_point() + xlab("Patient id") + ylab("Estradiol (nmol/L)")
figureSMG.a <- figureSMG.a + theme(text = element_text(size=15),
                                 axis.line = element_line(size = 1.25),
                                 axis.ticks = element_line(size = 2),
                                 axis.ticks.length=unit(.25, "cm"))

## ** survival curves
e.cox <- coxph(Surv(max(estra)-estra, observed) ~ strata(group),
               data = dt.estra, x = TRUE, y = TRUE)
e.pred <- predictCox(e.cox, keep.newdata = TRUE)

figureSMG.b <- autoplot(e.pred, group.by = "strata", type = "survival", plot = FALSE)
figureSMG.b <- figureSMG.b$plot + scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25),
                                                     labels = -c(0,0.25,0.5,0.75,1,1.25)+max(dt.estra$estra))
figureSMG.b <- figureSMG.b + ylab("Survival") + xlab("Estradiol (nmol/L)") + labs(colour="OC group", shape = "Left censored")
figureSMG.b <- figureSMG.b  + scale_shape_manual("Below the detection limit",
                                                 values = c(0,19),
                                                 labels = c("yes","no"))
figureSMG.b <- figureSMG.b + theme(text = element_text(size=15),
                                   axis.line = element_line(size = 1.25),
                                   axis.ticks = element_line(size = 2),
                                   axis.ticks.length=unit(.25, "cm"))

## ** assemble
figureSMG <- ggarrange(figureSMG.a, figureSMG.b + coord_cartesian(ylim = c(0,1)), ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE)

## * Export
ggsave(figureSMG, file = file.path(path.figures,"figureSMG.pdf"), width = 10, height = 6)



