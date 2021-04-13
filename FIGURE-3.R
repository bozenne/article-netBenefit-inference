path.figures <- "figures-article"

## * data + packages
library(ggplot2)
library(ggpubr)
library(flexsurv)
source("ILLUSTRATION-0-data-management.R")

## * BuyseTest
dt.prodige2 <- copy(dt.prodige)
dt.prodige2$effects <- dt.prodige2$toxicity
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(effects, threshold = 0.5, operator = "<0"),
                  data = dt.prodige2, method.inference = "u-statistic")
eSe.BT <- sensitivity(e.BT, threshold = list(OS=seq(0,12,1), effects = 1:3), band = FALSE, adj.p.value = FALSE)
test <- (eSe.BT$effects == 1) * (eSe.BT$OS == 2)
eSe.BT$color <- c("black","blue")[test+1]

figure3 <- autoplot(eSe.BT, plot = FALSE, col = NA, label = "Threshold for side") + geom_abline(slope=0, intercept = 0, color = "red") + xlab("Threshold for OS (in months)")
figure3 <- figure3 + geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), color = eSe.BT$color,size=1.5*test+0.75) + geom_point(color = eSe.BT$color, size = 2*test+2)
figure3 <- figure3 + guides(color = FALSE)
figure3 <- figure3 + theme(text = element_text(size=25),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"))
## * export
ggsave(figure3, filename = file.path(path.figures,"figure3.pdf"), height = 7, width = 14)

