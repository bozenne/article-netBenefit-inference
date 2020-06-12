### FCT-gg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  4 2020 (11:05) 
## Version: 
## Last-Updated: jun  6 2020 (15:06) 
##           By: Brice Ozenne
##     Update #: 217
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ggTiming
ggTiming <- function(data, file = NULL, plot = TRUE, txt.size = 20){
    require(ggpubr)

    Ufile <- unique(data$iFile)
    if(is.null(file)){
        file <- 1:length(Ufile)
    }else  if(any(file %in% 1:length(Ufile) == FALSE)){
        stop("unknown file selected \n")
    }
    if("threshold" %in% names(data) && length(unique(data$threshold))>1){
        stop("Cannot handle multple thresholds")
    }
    if("endpoint" %in% names(data) && length(unique(data$endpoint))>1){
        stop("Cannot handle multple endpoints")
    }

    dtW.gg <- data[method %in% c("Gehan","Peron") & Hprojection == 1 & iFile %in% Ufile[file],
                   .SD, .SDcols = c("n","method","timeAll","timeEstimate")]

    dtW.gg[, time := timeAll/timeEstimate]
    dtL1.gg <- dtW.gg[, .(n,method,time)]

    dtL2.gg <- melt(dtW.gg,
                    id.vars = c("n","method"),
                    measure.vars = c("timeAll","timeEstimate"),
                    value.name = "time",
                    variable.name = "type")

    dtL2.gg[, type := factor(type, levels = c("timeEstimate","timeAll"),
                             labels = c("estimate","estimate+se"))]
    dtL.gg <- rbind(cbind(dtL1.gg[,.(n,method)], type = "relative", dtL1.gg[,.(time)]),
                    cbind(dtL2.gg))
   
    dtL.gg[, n := factor(n, levels = sort(unique(n)))]
    dtL.gg[, scoring.rule := factor(method, levels = c("Gehan","Peron"),
                                    labels = c("Gehan's scoring rule","Peron's scoring rule"))]
    dtL.gg[, type := factor(type,
                            levels = c("estimate","estimate+se","relative"),
                            labels = c("estimate","estimate with s.e.","ratio"))]

    gg1 <- ggplot(data = dtL.gg[type == levels(type)[1]], mapping = aes(x = n, y = time))
    gg1 <- gg1 + geom_boxplot()
    gg1 <- gg1 + facet_grid(type~scoring.rule)
    gg1 <- gg1 + xlab("") + ylab("time (s)") + labs(color="")
    gg1 <- gg1 + theme(legend.position="bottom",
                       text = element_text(size=txt.size),
                       plot.margin = unit(c(0.1,0.1,0,0.35), "cm"))

    gg2 <- ggplot(data = dtL.gg[type == levels(type)[2]], mapping = aes(x = n, y = time))
    gg2 <- gg2 + geom_boxplot()
    gg2 <- gg2 + facet_grid(type~scoring.rule)
    gg2 <- gg2 + xlab("") + ylab("time (s)") + labs(color="")
    gg2 <- gg2 + theme(legend.position="bottom",
                       text = element_text(size=txt.size),
                       plot.margin = unit(c(-0.1,0.1,0.2,0), "cm"))

    gg3 <- ggplot(data = dtL.gg[type == levels(type)[3]], mapping = aes(x = n, y = time))
    gg3 <- gg3 + geom_hline(yintercept = 1, color = "red", size = 1.5)
    gg3 <- gg3 + geom_boxplot()
    gg3 <- gg3 + facet_grid(type~scoring.rule)
    gg3 <- gg3 + xlab("sample size in each group") + ylab("relative time") + labs(color="")
    gg3 <- gg3 + theme(legend.position="bottom",
                       text = element_text(size=txt.size),
                       plot.margin = unit(c(-0.1,0.1,0.2,0.5), "cm"))

    gg <- ggarrange(gg1, gg2, gg3, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")

    if(plot){
        print(gg)
    }
    
    return(invisible(list(plot = gg,
                          data = dtL.gg)))
}

## * ggBias
ggBias <- function(data, file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    Ufile <- unique(data$iFile)
    if(is.null(file)){
        file <- length(Ufile)
    }else  if(any(file %in% 1:length(Ufile) == FALSE)){
        stop("unknown file selected \n")
    }
    if("threshold" %in% names(data) && length(unique(data$threshold))>1){
        stop("Cannot handle multple thresholds")
    }
    if("endpoint" %in% names(data) && length(unique(data$endpoint))>1){
        stop("Cannot handle multple endpoints")
    }

    dtW.gg <- data[Hprojection == 1 & iFile %in% Ufile[file],.(n,method,estimate)]
    dtW.gg[, expected := mean(.SD[n==max(n) & method == "GS",estimate])]

    dtW.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
    dtW.gg[, n := factor(n, levels = sort(unique(n)))]    
    dtW.gg[, sample.size := factor(paste0("sample size = ",n), paste0("sample size = ",unique(n)))]

    gg <- ggplot(dtW.gg, aes(y = estimate - expected))
    gg <- gg + facet_wrap(~sample.size)
    gg <- gg + geom_hline(yintercept = 0, color = "red", size = 1.5)
    gg <- gg + geom_boxplot(aes(color = scoring.rule))#, outlier.shape = NA)
    gg <- gg + xlab("sample size") + ylab("bias") + labs(color="")
    gg <- gg + theme(legend.position="bottom",
                     text = element_text(size=txt.size),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

    if(plot){
        print(gg)
    }
    
    return(invisible(list(plot = gg,
                          data = dtW.gg)))
}

## * ggSe
ggSe <- function(data, file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    Ufile <- unique(data$iFile)
    if(is.null(file)){
        file <- 1:length(Ufile)
    }else  if(any(file %in% 1:length(Ufile) == FALSE)){
        stop("unknown file selected \n")
    }
    if("threshold" %in% names(data) && length(unique(data$threshold))>1){
        stop("Cannot handle multple thresholds")
    }
    if("endpoint" %in% names(data) && length(unique(data$endpoint))>1){
        stop("Cannot handle multple endpoints")
    }

    dtW.gg <- data[iFile %in% Ufile[file],.(n,method,Hprojection,se,estimate)]
    dtS.gg <- dtW.gg[,.(model = mean(se), empirical = sd(estimate)),by = c("n","method","Hprojection")]
    dtS.gg[, n := factor(n, levels = sort(unique(n)))]
    dtS.gg[, method2 := paste0(method," H",Hprojection)]
    dtS.gg[, Hprojection := factor(Hprojection, levels = 1:2, labels = c("1st order H-projection","2nd order H-projection"))]
    dtS.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
    
    gg <- ggplot(dtS.gg, aes(x = empirical, y = model, group = Hprojection, color = n))
    gg <- gg + geom_abline(slope = 1)
    gg <- gg + geom_point(size = 2) + geom_line(size = 1.25)
    gg <- gg + facet_grid(Hprojection~scoring.rule)
    gg <- gg + xlab("empirical standard error") + ylab("average estimated standard error") 
    gg <- gg + labs(colour = "sample size (per group)") + theme(legend.position="bottom", text = element_text(size=txt.size))

    if(plot){
        print(gg)
    }
    
    return(invisible(list(plot = gg,
                          data = dtS.gg)))

}

## * ggCoverage
ggCoverage <- function(data, file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    require(ggthemes)
    
    Ufile <- unique(data$iFile)
    if(is.null(file)){
        file <- 1:length(Ufile)
    }else  if(any(file %in% 1:length(Ufile) == FALSE)){
        stop("unknown file selected \n")
    }

    keep.cols <- c("n", "method", "Hprojection", "threshold", "endpoint", "estimate", "lower.ci", "upper.ci")
    keep.cols <- keep.cols[keep.cols %in% names(data)]

    test.by <- sapply(c("threshold","endpoint","Hprojection"), function(iBy){
        length(unique(data[[iBy]]))>1
    })
    if(sum(test.by)>1){
        stop("Only one of \"threshold\", \"endpoint\", and \"Hprojection\" can vary \n")
    }
    if(sum(test.by)>0){
        by <- names(test.by)[which(test.by)]
    }else{
        by <- NULL
    }
    
    dtW.gg <- data[iFile %in% Ufile[file],.SD,.SDcols = keep.cols]
    if(is.null(expected)){
        dtW.gg[, expected := mean(.SD[n==max(n) & method == "GS",estimate]), by = by]
    }else{
        Uby <- unique(dtW.gg[[by]])
        n.by <- length(Uby)
        for(iby in 1:n.by){
            dtW.gg[dtW.gg[[by]] == Uby[iby], expected := expected[iby]]
        }
    }
    dtW.gg[, coverage := (lower.ci < expected)*(upper.ci > expected)]
    dtS.gg <- dtW.gg[, .(rep = .N, coverage =  mean(coverage)), by = c("n","method","Hprojection",by)]
    
    dtS.gg[, n := factor(n, levels = sort(unique(n)))]
    dtS.gg[, method2 := paste0(method," H",Hprojection)]
    dtS.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
    if(identical(by,"Hprojection")){
        dtS.gg[, Hprojection := factor(Hprojection, levels = 1:2, labels = c("1st order H-projection","2nd order H-projection"))]
    }else if(identical(by,"threshold")){
        dtS.gg[, threshold := paste0("\u03C4=", threshold)]
    }else if(identical(by,"endpoint")){

    }
    
    gg <- ggplot(dtS.gg, aes(x = n, y = coverage, group = scoring.rule, color = scoring.rule))
    gg <- gg + geom_hline(yintercept = 0.95, color = "red", size = 1.5)
    gg <- gg + geom_point(size = 2) + geom_line(size = 1.5)
    if(length(by)>0){
        gg <- gg + facet_grid(as.formula(paste0("~",by)))
    }
    gg <- gg + theme(legend.position="bottom", text = element_text(size=txt.size),
                     axis.text.x = element_text(angle = 90, hjust = 1))
    gg <- gg + xlab("sample size in each group") + ylab("coverage") + labs(colour="")
    gg <- gg + scale_colour_colorblind()

    if(plot){
        print(gg)
    }
    
    return(invisible(list(plot = gg,
                          data = dtS.gg)))

}

## * createTable
createTable <- function(data, file = NULL, expected = NULL,
                        print = TRUE, digits = NULL,
                        label = "", caption = "", trace = TRUE){
    require(xtable)

    ## ** select data
    Ufile <- unique(data$iFile)
    if(is.null(file)){
        file <- 1:length(Ufile)
    }else  if(any(file %in% 1:length(Ufile) == FALSE)){
        stop("unknown file selected \n")
    }
    dt.table <- data[iFile %in% Ufile[file],.(n,method, Hprojection, threshold,
                                              estimate, se, lower.ci, upper.ci,
                                              timeEstimate, timeAll)]
    if(is.null(expected)){
        dt.table[, expected := mean(.SD[n==max(n) & method == "GS",estimate]), by = "threshold"]
    }else{
        Uthreshold <- unique(dt.table$threshold)
        n.threshold <- length(Uthreshold)
        for(iTh in 1:n.threshold)
            dt.table[threshold == Uthreshold[iTh], expected := expected[iTh]]
    }


    ## ** summarize
    dt.table[, coverage := (lower.ci < expected)*(upper.ci>expected)]
    byVar <- c("n"[length(unique(dt.table$n))>1],
               "Hprojection"[length(unique(dt.table$Hprojection))>1],
               "threshold"[length(unique(dt.table$threshold))>1],
               "method"[length(unique(dt.table$method))>1]
               )
    if("Hprojection" %in% byVar){stop("Cannot handle several Hprojections")}
    dtS.table <- dt.table[, .(rep = .N,
                              bias = mean(estimate - expected),
                              empirical= sd(estimate),
                              estimated= mean(se),
                              coverage =  mean(coverage)),
                          by = byVar]
    
    if(trace){cat("Number of repetitions: ",paste(unique(dtS.table$rep),collapse = " "),"\n")}
    dtS.table[, rep := NULL]

    ## ** create table
    setkeyv(dtS.table,c("n","threshold"))
    dtS.table$threshold[duplicated(interaction(dtS.table$n,dtS.table$threshold))] <- ""
    dtS.table$n <- paste0("\\(n=m=\\) ",dtS.table$n)
    dtS.table$n[duplicated(dtS.table$n)] <- ""
    
    setnames(dtS.table, old = "n", new = "sample size")

    setnames(dtS.table, old = "threshold", new = "\\(\\tau\\)")
        dtS.table[, method := factor(method, levels = c("GS","Gehan","Peron"),
                                 labels = c("Full data","Gehan", "Peron"))]
    setnames(dtS.table, old = "method", new = "scoring rule")
    setnames(dtS.table, old = "empirical", new = "empirical \\(\\sigma_{\\hat{\\Delta},\\hat{\\Delta}}\\)")
    setnames(dtS.table, old = "estimated", new = "estimated \\(\\sigma_{\\hat{\\Delta},\\hat{\\Delta}}\\)")

    out <- xtable(dtS.table, type = "latex", 
                  label = label,
                  caption = caption,
                  digits = digits)


    n.method <- length(unique(dt.table$method))
    addtorow <- list()
    addtorow$pos <- as.list(seq(from = n.method, to = NROW(dtS.table), by=n.method))
    addtorow$command <- rep(c(rep("[2mm]", times = n.method-1),"[4mm]"),
                            times = length(addtorow$pos)/n.method)

    addtorow$pos <- addtorow$pos[-length(addtorow$pos)]
    addtorow$command <- addtorow$command[-length(addtorow$command)]
    
    if(print){ print(out, add.to.row = addtorow, include.rownames=FALSE, include.colnames = TRUE,
                     sanitize.colnames.function = identity,
                     sanitize.text.function = identity,
                     table.placement = "!h")}
    return(invisible(out))
}

######################################################################
### FCT-gg.R ends here
