path.figures <- "figures-article"

## * data + packages
source("ILLUSTRATION-0-data-management.R")

## * BuyseTest
e.BT <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2) + cont(toxicity, threshold = 1, operator = "<0"),
                  data = dt.prodige, method.inference = "u-statistic")

e.BT_boot <- BuyseTest(bras ~ tte(OS, status = etat, threshold = 2)  + cont(toxicity, threshold = 1, operator = "<0"),
                       data = dt.prodige, method.inference = "bootstrap", n.resampling = 1e4, seed = 10, cpus = 10,
                       trace = 1)

table3 <- xtable(apply(summary(e.BT, print = FALSE, digit = c(2,4,2))$table.print[,-c(2:3,12),],2,as.character))
print(table3, include.rownames = FALSE)
## \begin{table}[ht]
## \centering
## \begin{tabular}{lllllllll}
##   \hline
## endpoint & favorable(\%) & unfavorable(\%) & neutral(\%) & uninf(\%) & delta & Delta & CI [2.5\% ; 97.5\%] & p.value \\ 
##   \hline
## OS & 56.88 & 26.49 & 16.61 & 0.02 &  0.3039 & 0.3039 & [0.1819;0.4167] & 2.2e-06 \\ 
##   toxicity &  4.41 &  6.53 &  5.69 & 0.00 & -0.0212 & 0.2827 & [0.1582;0.3983] & 1.4e-05 \\ 
##    \hline
## \end{tabular}
## \end{table}
