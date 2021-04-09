path.figures <- "figures-article"

## * data + packages
library(ggplot2)
library(ggpubr)
library(flexsurv)
source("ILLUSTRATION-0-data-management.R")

## * BuyseTest
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 0.5, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic")
eSe.BT <- sensitivity(e.BT, threshold = list(OS=seq(0,12,1), toxicity = 1:3), band = TRUE, adj.p.value = FALSE)

figure4 <- autoplot(eSe.BT, plot = FALSE, col = NA) + geom_abline(slope=0, intercept = 0, color = "red") + xlab("Threshold for OS (in months)")
figure4 <- figure4 + theme(text = element_text(size=25),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"))
## * export
ggsave(figure4, filename = file.path(path.figures,"figure4.pdf"), height = 7, width = 14)

