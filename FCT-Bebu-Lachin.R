### FCT-Bebu-Lachin.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  7 2020 (09:24) 
## Version: 
## Last-Updated: jul  7 2020 (11:13) 
##           By: Brice Ozenne
##     Update #: 46
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
UstatBL <- function(formula, data, method){
    require(data.table)

    ## ** extract information
    var.group <- all.vars(update(formula,".~0"))
    var.endpoint <- all.vars(update(formula,"0~."))
    Ugroup <- unique(data[[var.group]])

    ## ** check arguments
    method <- match.arg(method, choices = c("original","corrected","alternative"))
    if(length(var.endpoint)>1){stop("only valid for a single endpoint \n")}
    if(length(Ugroup)!=2){stop("the group variable must take exactly 2 different values \n")}
    if(!is.numeric(data[[var.endpoint]])){stop("the endpoint variable must be numeric \n")}

    ## ** normalize arguments
    X <- data[[var.endpoint]][data[[var.group]]==Ugroup[2]]
    m <- length(X)
    Y <- data[[var.endpoint]][data[[var.group]]==Ugroup[1]]
    n <- length(Y)
    N <- m + n

    ## ** sufficient statistics
    grid <- cbind(expand.grid(i = 1:m, j = 1:n),expand.grid(Xi = X,Yj = Y))
    grid$Sij1 <- apply(grid,1, function(x){x["Xi"]>x["Yj"]})
    grid$Sij2 <- apply(grid,1, function(x){x["Xi"]<x["Yj"]})
    dt.grid <- as.data.table(grid)

    Si <- dt.grid[,.(Si1 = sum(Sij1), Si2 = sum(Sij2)),by = "i"]
    Sj <- dt.grid[,.(Sj1 = sum(Sij1), Sj2 = sum(Sij2)),by = "j"]

    ## ** variance calculation
    pXY1 <- sum(dt.grid$Sij1) / (m*n) ## P[X1>Y1] = 1/mn \sum_ij S1_ij
    pXY2 <- sum(dt.grid$Sij2) / (m*n) ## P[X1<Y1] = 1/mn \sum_ij S2_ij

    if(method == "original"){
        pXY1XYp1 <- sum(Si$Si1*(Si$Si1-1)) / (m*n*(n-1)) ## P[X1>Y1,X1>Y1'] = 1/(mn(n-1)) \sum_ij S1_i. (S1_i. - 1)
        pXY1XpY1 <- sum(Sj$Sj1*(Sj$Sj1-1)) / (n*m*(m-1)) ## P[X1>Y1,X1'>Y1] = 1/(mn(n-1)) \sum_ij S1_.j (S1_.j - 1)
        pXY2XYp2 <- sum(Si$Si2*(Si$Si2-1)) / (m*n*(n-1)) ## P[X1<Y1,X1<Y1'] = 1/(mn(n-1)) \sum_ij S2_i. (S2_i. - 1)
        pXY2XpY2 <- sum(Sj$Sj2*(Sj$Sj2-1)) / (n*m*(m-1)) ## P[X1<Y1,X1'<Y1] = 1/(mn(n-1)) \sum_ij S2_.j (S2_.j - 1)
        pXY1XYp2 <- sum(Si$Si1*Si$Si2) / (m*n*(n-1)) ## P[X1>Y1,X1<Y1'] = 1/(mn(n-1)) \sum_ij S1_i. S2_i.
        pXY1XpY2 <- sum(Sj$Sj1*Sj$Sj2) / (n*m*(m-1)) ## P[X1>Y1,X1'<Y1] = 1/(mn(n-1)) \sum_ij S1_.j S2_.j
    }else{
        gridI <- cbind(expand.grid(i1 = 1:m, i2 = 1:m, j = 1:n), expand.grid(Xi1 = X, Xi2 = X, Yj = Y))
        gridJ <- cbind(expand.grid(i = 1:m, j1 = 1:n, j2 = 1:n), expand.grid(Xi = X, Yj1 = Y, Yj2 = Y))
        if(method == "alternative"){
            gridI <- gridI[gridI$i1!=gridI$i2,]
            gridJ <- gridJ[gridJ$j1!=gridJ$j2,]
        }
        pXY1XYp1 <- mean(apply(gridJ,1, function(x){(x["Xi"]>x["Yj1"])*(x["Xi"]>x["Yj2"])}))
        pXY1XpY1 <- mean(apply(gridI,1, function(x){(x["Xi1"]>x["Yj"])*(x["Xi2"]>x["Yj"])}))
        pXY2XYp2 <- mean(apply(gridJ,1, function(x){(x["Xi"]<x["Yj1"])*(x["Xi"]<x["Yj2"])}))
        pXY2XpY2 <- mean(apply(gridI,1, function(x){(x["Xi1"]<x["Yj"])*(x["Xi2"]<x["Yj"])}))
        pXY1XYp2 <- mean(apply(gridJ,1, function(x){(x["Xi"]>x["Yj1"])*(x["Xi"]<x["Yj2"])}))
        pXY1XpY2 <- mean(apply(gridI,1, function(x){(x["Xi1"]>x["Yj"])*(x["Xi2"]<x["Yj"])}))
    }
    
    xi_11_10 <- pXY1XYp1 - pXY1^2 ## \Cov(X>Y,X>Y') = \Esp[(X>Y)(X>Y')] - \Esp[X>Y]\Esp[X>Y']
    xi_22_10 <- pXY2XYp2 - pXY2^2 ## \Cov(X<Y,X<Y') = \Esp[(X<Y)(X<Y')] - \Esp[X<Y]\Esp[X<Y']
    xi_12_10 <- pXY1XYp2 - pXY1*pXY2 ## \Cov(X>Y,X<Y') = \Esp[(X>Y)(X<Y')] - \Esp[X>Y]\Esp[X<Y']
    xi_21_10 <- pXY1XYp2 - pXY2*pXY1 ## \Cov(X<Y,X>Y') = \Esp[(X<Y)(X>Y')] - \Esp[X<Y]\Esp[X>Y']

    xi_11_01 <- pXY1XpY1 - pXY1^2 ## \Cov(X>Y,X'>Y) = \Esp[(X>Y)(X'>Y)] - \Esp[X>Y]\Esp[X'>Y]
    xi_22_01 <- pXY2XpY2 - pXY2^2 ## \Cov(X<Y,X'<Y) = \Esp[(X<Y)(X'<Y)] - \Esp[X<Y]\Esp[X'<Y]
    xi_12_01 <- pXY1XpY2 - pXY1*pXY2 ## \Cov(X>Y,X'<Y) = \Esp[(X>Y)(X'<Y)] - \Esp[X>Y]\Esp[X'<Y]
    xi_21_01 <- pXY1XpY2 - pXY2*pXY1 ## \Cov(X<Y,X'>Y) = \Esp[(X<Y)(X'>Y)] - \Esp[X<Y]\Esp[X'>Y]

    sigma_11 <- (N/m)*xi_11_10 + (N/n)*xi_11_01
    sigma_22 <- (N/m)*xi_22_10 + (N/n)*xi_22_01
    sigma_12 <- (N/m)*xi_12_10 + (N/n)*xi_12_01
    sigma_21 <- (N/m)*xi_21_10 + (N/n)*xi_21_01

    ## ** export
    return(list(estimate = c(favorable = pXY1, unfavorable = pXY2, netBenefit = pXY1 - pXY2),
                sd = c(favorable = sqrt(sigma_11/N), unfavorable = sqrt(sigma_22/N), netBenefit = sqrt((sigma_11 + sigma_22 - sigma_12 - sigma_21)/N)),
                variance = matrix(c(sigma_11,sigma_12,sigma_21,sigma_22)/N, nrow = 2, ncol = 2,
                                  dimnames = list(c("favorable","unfavorable"),
                                                  c("favorable","unfavorable")))
                ))
}

######################################################################
### FCT-Bebu-Lachin.R ends here
