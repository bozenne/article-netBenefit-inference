## * data + packages
library(ggplot2)
source("ILLUSTRATION-data-management.R")

## * fit
## ** parametric
AFT0 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1 <- flexsurvreg(Surv(OS, etat) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0 <- exp(coef(AFT0))
param1 <- exp(coef(AFT1))

AFT0.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Gemcitabine",], dist = "Weibull")
AFT1.cens <- flexsurvreg(Surv(OS, etat==0) ~ 1, data = dt.prodige[dt.prodige$bras == "Folfirinox",], dist = "Weibull")
param0.cens <- exp(coef(AFT0.cens))
param1.cens <- exp(coef(AFT1.cens))

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

## ** non-parametric
e.coxph <- coxph(Surv(OS,etat)~strata(bras), data = dt.prodige, x = TRUE, y = TRUE)
pred.coxph <- predictCox(e.coxph, type = "survival", keep.newdata = TRUE)

## * display
figure2 <- autoplot(pred.coxph,plot = FALSE)$plot

figure2 <- figure2 + geom_line(data = dt, mapping = aes(x = time, y = surv, group = strata, color = type), size = 0.8)
figure2$layers <- figure2$layers[c(3,1,2)]
figure2 <- figure2 + scale_colour_manual(name = "",
                                           values = c("Weibull model" = "black",
                                                      "Gemcitabine" = "darkorange",
                                                      "Folfirinox" = "darkblue"
                                                      ),
                                           labels = c("Weibull model" = "Weibull model",
                                                      "Gemcitabine" = "arm Folfirinox",
                                                      "Folfirinox" = "arm Gemcitabine"
                                                      ))
figure2 <- figure2 + scale_shape_manual(name = "",
                                          breaks = c(0,1),
                                          values = c(3,18),
                                          labels = c("censoring","death"))
figure2 <- figure2 + xlab("Months") + ylab("Probability of survival")
figure2 <- figure2 + theme(text = element_text(size=25),
                           axis.line = element_line(size = 1.25),
                           axis.ticks = element_line(size = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.box.margin = margin(c(-40,0,0,0)),
                           legend.position="bottom",
                           legend.direction = "vertical",
                           legend.background = element_rect(fill = NA))
figure2

if(FALSE){
    ggsave(figure2, filename = file.path("./figures/ILLUSTRATION-folfirinox-surv.pdf"), height = 7, width = 7.5)
}
