path.figures <- "figures-article"

## * data + packages
library(ggplot2)
library(ggpubr)
library(flexsurv)
source("ILLUSTRATION-0-data-management.R")

## * prepare
## ** efficacy: fit parametric models
AFT0 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0 <- exp(coef(AFT0))
param1 <- exp(coef(AFT1))

X <- seq(0,50, by = 0.1)
dt <- rbind(data.table(strata = "Gemcitabine",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param0["scale"], shape = param0["shape"])
                       ),
            data.table(strata = "Folfirinox",
                       time = X,
                       type = "Weibull model",
                       surv = 1-pweibull(q = X, scale = param1["scale"], shape = param1["shape"])
                       )
            )

## ** efficacy: fit non-parametric model
e.coxph <- coxph(Surv(OS,etat)~strata(bras), data = dt.prodige, x = TRUE, y = TRUE)
pred.coxph <- predictCox(e.coxph, type = "survival", keep.newdata = TRUE)

## plot(prodlim(Hist(OS,etat)~bras, data = dt.prodige))

## ** toxicity
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

## ** efficacy: survival curves
figure2a <- autoplot(pred.coxph, group.by = "strata", plot = FALSE)$plot

figure2a <- figure2a + geom_line(data = dt, mapping = aes(x = time, y = surv, group = strata, color = type), size = 0.8)
figure2a$layers <- figure2a$layers[c(3,1,2)]
figure2a <- figure2a + scale_colour_manual(name = "",
                                           values = c("Folfirinox" = "darkblue",
                                                      "Gemcitabine" = "darkorange",
                                                      "Weibull model" = "black"                                                      
                                                      ),
                                           labels = c("Folfirinox" = "arm Folfirinox",
                                                      "Gemcitabine" = "arm Gemcitabine",
                                                      "Weibull model" = "Weibull model"
                                                      ))
figure2a <- figure2a + scale_shape_manual(name = "",
                                          breaks = c(0,1),
                                          values = c(3,18),
                                          labels = c("censoring","death"))
figure2a <- figure2a + xlab("Months") + ylab("Probability of survival")
figure2a <- figure2a + theme(text = element_text(size=25),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.box.margin = margin(c(-40,0,0,0)),
                           legend.position="bottom",
                           legend.direction = "vertical",
                           legend.background = element_rect(fill = NA))
## figure2a

## ** toxiciy: barplot
figure2b <- ggplot(dt.toxicity, aes(fill=toxicity, y = pc, x = bras))
figure2b <- figure2b + geom_col(position = position_stack(reverse = TRUE))
figure2b <- figure2b + scale_fill_manual("Worst adverse event", values = cc)
figure2b <- figure2b + xlab("") + ylab("Relative frequency") 
figure2b <- figure2b + scale_y_continuous(labels = scales::percent) 
figure2b <- figure2b + theme(text = element_text(size=25),
                           axis.text.x = element_text(colour = c("darkblue","orange")),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.box.margin = margin(c(-35,0,0,0)),
                           legend.position="bottom",
                           legend.direction = "vertical")
figure2b <- figure2b + guides(fill=guide_legend(nrow=1,byrow=TRUE))
## figure2b

## ** assemble
figure2 <- ggarrange(figure2a, figure2b)

## * export
ggsave(figure2, filename = file.path(path.figures,"figure2.pdf"), height = 7, width = 14)


