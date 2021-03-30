## * data + packages
library(ggplot2)
source("ILLUSTRATION-data-management.R")

## * compute percentage
dt.toxicity <- dt.prodige[,.N, by = c("toxicity","bras")]
dt.toxicity[, bras := paste0("arm \n",bras)]
dt.toxicity[, Ntot := sum(N), by = "bras"]
dt.toxicity[, pc := N/Ntot]
setkeyv(dt.toxicity, c("bras","toxicity"))
dt.toxicity[, toxicity := factor(toxicity, levels = 0:5)]

cc <- scales::seq_gradient_pal(rgb(r=0,g=0.9,b=0),
                               rgb(r=0.9,g=0,b=0),
                               "Lab")(seq(0,1,length.out=6))

## * display
figure3 <- ggplot(dt.toxicity, aes(fill=toxicity, y = pc, x = bras))
figure3 <- figure3 + geom_col(position = position_stack(reverse = TRUE))
figure3 <- figure3 + scale_fill_manual("Worst adverse event", values = cc)
figure3 <- figure3 + xlab("") + ylab("Relative frequency") 
figure3 <- figure3 + scale_y_continuous(labels = scales::percent) 
figure3 <- figure3 + theme(text = element_text(size=25),
                           axis.text.x = element_text(colour = c("orange","darkblue")),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.box.margin = margin(c(-35,0,0,0)),
                           legend.position="bottom",
                           legend.direction = "vertical")
figure3 <- figure3 + guides(fill=guide_legend(nrow=1,byrow=TRUE))
figure3
## ggsave(figure3, filename = "figures/ILLUSTRATION-folfirinox-tox.pdf", width = 6, height = 7)
