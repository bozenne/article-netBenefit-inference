## * packages
library(BuyseTest)
library(survival)
library(prodlim)
library(flexsurv)
library(riskRegression)
library(ggplot2)
library(ggthemes)
library(data.table)

## * data management
## ** load
dt.prodige <- as.data.table(read.csv(file.path("data","PRODIGE.csv"), sep = ";", header = TRUE))

## ** process
dt.prodige[, d_dn2 := as.Date(d_dn, "%d/%m/%Y")]
dt.prodige[, randodt2 := as.Date(randodt, "%d/%m/%Y")]
dt.prodige[, bras := factor(ifelse(bras==2,"Gemcitabine","Folfirinox"), c("Gemcitabine","Folfirinox"))]
dt.prodige[, OS := as.numeric(difftime(d_dn2,randodt2,units="days")/30.44)]

M.tox <- cbind(inf = dt.prodige[,pmax(mx_inf1, mx_inf2, mx_inf3, mx_inf4, 0, na.rm = TRUE)],
               car = dt.prodige[,pmax(mx_car1, mx_car2, mx_car3, mx_car4, mx_car5, mx_car6, mx_car7, 0, na.rm = TRUE)],
               aon = dt.prodige[,pmax(mx_aon1, mx_aon2, mx_aon3, mx_aon4, mx_aon5, mx_aon6, mx_aon7, 0, na.rm = TRUE)],
               pul = dt.prodige[,pmax(mx_pul1, mx_pul2, mx_pul3, mx_pul4, 0, na.rm = TRUE)],
               aut = dt.prodige[,pmax(mx_aut1, mx_aut2, mx_aut3, mx_aut4, 0, na.rm = TRUE)],
               ren = dt.prodige[,pmax(mx_ren1, mx_ren2, mx_ren3, mx_ren4, 0, na.rm = TRUE)],
               gas = dt.prodige[,pmax(mx_gas1, mx_gas2, mx_gas3, mx_gas4, mx_gas5, mx_gas6, mx_gas7, 0, na.rm = TRUE)],
               der = dt.prodige[,pmax(mx_der1, mx_der2, mx_der3, mx_der4, mx_der5, mx_der6, mx_der7, 0, na.rm = TRUE)],
               tog = dt.prodige[,pmax(mx_tog1, mx_tog2, mx_tog3, mx_tog4, 0, na.rm = TRUE)],
               uro = dt.prodige[,pmax(mx_uro1, mx_uro2, mx_uro3, mx_uro4, 0, na.rm = TRUE)]
               )
dt.prodige[,toxicity := apply(M.tox,1, max)]
dt.prodige[, strata := interaction(oms_r,loca_r)]
  
rm.cols <- setdiff(names(dt.prodige),c("OS","bras","toxicity","etat","strata"))
dt.prodige[,c(rm.cols) := NULL]
