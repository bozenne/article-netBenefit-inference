### FCT-gg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  4 2020 (11:05) 
## Version: 
## Last-Updated: Apr 12 2021 (10:40) 
##           By: Brice Ozenne
##     Update #: 299
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ggTiming
ggTiming <- function(data, type.data = "raw", file = NULL, plot = TRUE, txt.size = 20){
    require(ggpubr)
    type.data <- match.arg(type.data, c("raw","processed"))

    if(type.data == "raw"){
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
    }else if(type.data == "processed"){
        dtL.gg <- data
    }

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

    if(file.exists("tempo.pdf")){
        stop("Cannot run the function as a file name tempo.pdf is in the working directory. \n")
    }
    pdf("tempo.pdf")
    gg <- try(ggarrange(plotlist = list(gg1, gg2, gg3), ncol=1, nrow=3, common.legend = TRUE, legend="bottom"), silent = TRUE)
    dev.off()
    file.remove("tempo.pdf")
    
    if(plot){
        print(gg)
    }
    
    return(invisible(list(plot = gg,
                          data = dtL.gg)))
}

## * ggBias
ggBias <- function(data, type.data = "raw", file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    type.data <- match.arg(type.data, c("raw","processed"))

    if(type.data == "raw"){
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
        dtW.gg <- data[Hprojection == 1 & iFile %in% Ufile[file],.(n,method,estimate)]
        if(is.numeric(dtW.gg$n)){
            max.n <- max(dtW.gg$n)
        }else if(is.factor(dtW.gg$n)){
            max.n <- tail(levels(dtW.gg$n),1)
        }
        dtW.gg[, expected := mean(.SD[n==max.n & method == "GS",estimate])]
        dtW.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
        dtW.gg[, n := factor(n, levels = sort(unique(n)))]    
        dtW.gg[, sample.size := factor(paste0("sample size = ",n), paste0("sample size = ",unique(n)))]
    }else if(type.data == "processed"){
        dtW.gg <- data
    }

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
ggSe <- function(data, type.data = "raw", file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    type.data <- match.arg(type.data, c("raw","processed"))

    if(type.data == "raw"){
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
        dtS.gg <- dtW.gg[,.(rep = .N, model = mean(se), empirical = sd(estimate)),by = c("n","method","Hprojection")]
        dtS.gg[, n := factor(n, levels = sort(unique(n)))]
        dtS.gg[, method2 := paste0(method," H",Hprojection)]
        dtS.gg[, Hprojection := factor(Hprojection, levels = 1:2, labels = c("1st order H-projection","2nd order H-projection"))]
        dtS.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
    }else if(type.data == "processed"){
        dtS.gg <- data
    }
    
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
ggCoverage <- function(data, by, type.data = "raw", file = NULL, plot = TRUE, expected = NULL, txt.size = 20){
    require(ggthemes)
    type.data <- match.arg(type.data, c("raw","processed"))
    
    if(type.data == "raw"){
        Ufile <- unique(data$iFile)
        if(is.null(file)){
            file <- 1:length(Ufile)
        }else  if(any(file %in% 1:length(Ufile) == FALSE)){
            stop("unknown file selected \n")
        }

        keep.cols <- c("n", "method", "Hprojection", "threshold", "endpoint", "estimate", "lower.ci", "upper.ci")
        keep.cols <- keep.cols[keep.cols %in% names(data)]

        dtW.gg <- data[iFile %in% Ufile[file],.SD,.SDcols = keep.cols]
        if(is.null(expected)){
            if(is.numeric(dtW.gg$n)){
                max.n <- max(dtW.gg$n)
            }else if(is.factor(dtW.gg$n)){
                max.n <- tail(levels(dtW.gg$n),1)
            }
            dtW.gg[, expected := mean(.SD[n==max.n & method == "GS",estimate]), by = by]
        }else{
            Uby <- unique(dtW.gg[[by]])
            n.by <- length(Uby)
            for(iby in 1:n.by){
                dtW.gg[dtW.gg[[by]] == Uby[iby], expected := expected[iby]]
            }
        }
        dtW.gg[, coverage := (lower.ci < expected)*(upper.ci > expected)]
        dtS.gg <- dtW.gg[, .(rep = .N, coverage =  mean(coverage, na.rm=TRUE)), by = c("n","method","Hprojection",by)]
    
        dtS.gg[, n := factor(n, levels = sort(unique(n)))]
        dtS.gg[, method2 := paste0(method," H",Hprojection)]
        dtS.gg[, scoring.rule := factor(method, levels = c("GS","Gehan","Peron"), labels = c("no censoring","Gehan's scoring rule","Peron's scoring rule"))]
        if(identical(by,"Hprojection")){
            dtS.gg[, Hprojection := factor(Hprojection, levels = 1:2, labels = c("1st order H-projection","2nd order H-projection"))]
        }else if(identical(by,"threshold")){
            dtS.gg[, threshold := paste0("\u03C4=", threshold)]
        }else if(identical(by,"endpoint")){

        }
    }else if(type.data == "processed"){
        dtS.gg <- data
    }

    test.by <- sapply(c("threshold","endpoint","Hprojection"), function(iBy){
        length(unique(dtS.gg[[iBy]]))>1
    })
    if(sum(test.by)>1){
        stop("Only one of \"threshold\", \"endpoint\", and \"Hprojection\" can vary \n")
    }
    if(sum(test.by)>0){
        by <- names(test.by)[which(test.by)]
    }else{
        by <- NULL
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
createTable <- function(data, by, type.data = "raw", file = NULL, expected = NULL,
                        print = TRUE, digits = NULL,
                        label = "", caption = "", trace = TRUE){
    require(xtable)
    type.data <- match.arg(type.data, c("raw","processed"))

    if(type.data == "raw"){
        ## ** select data
        Ufile <- unique(data$iFile)
        if(is.null(file)){
            file <- 1:length(Ufile)
        }else  if(any(file %in% 1:length(Ufile) == FALSE)){
            stop("unknown file selected \n")
        }
        keep.col <- c("n", "method", "Hprojection", by, "estimate", "se", "lower.ci", "upper.ci", "timeEstimate", "timeAll")
        dt.table <- data[iFile %in% Ufile[file], .SD, .SDcols = keep.col]
        if(is.null(expected)){
            if(is.numeric(dt.table$n)){
                max.n <- max(dt.table$n)
            }else if(is.factor(dt.table$n)){
                max.n <- tail(levels(droplevels(dt.table$n)),1)
            }
            dt.table[, expected := mean(.SD[n==max.n & method == "GS",estimate]), by = by]
        }else{
            Uby <- unique(dt.table[[by]])
            n.by <- length(Uby)
            for(iBy in 1:n.by)
                dt.table[dt.table[[by]] == Uby[iBy], expected := expected[iTh]]
        }

        ## ** summarize
        dt.table[, coverage := (lower.ci < expected)*(upper.ci>expected)]
        byVar <- c("n"[length(unique(dt.table$n))>1],
                   "Hprojection"[length(unique(dt.table$Hprojection))>1],
                   by[length(unique(dt.table[[by]]))>1],
                   "method"[length(unique(dt.table$method))>1]
                   )
        if("Hprojection" %in% byVar){stop("Cannot handle several Hprojections")}
        dtS.table <- dt.table[, .(rep = .N,
                                  seNA = sum(is.na(se)),
                                  bias = mean(estimate - expected),
                                  empirical= sd(estimate),
                                  estimated= mean(se, na.rm = TRUE),
                                  coverage =  mean(coverage, na.rm = TRUE)),
                              by = byVar]
    
    }else if(type.data == "processed"){
        dtS.table <- data.table::copy(data)
    }
    n.method <- length(unique(dtS.table$method))

    ## ** create table
    dtSS.table <- copy(dtS.table)
    if(trace){cat("Number of repetitions: ",paste(unique(dtSS.table$rep),collapse = " "),"\n")}
    dtSS.table[, rep := NULL]
    
    setkeyv(dtSS.table,c("n",by))
    dtSS.table[, seNA := NULL]
    
    dtSS.table[[by]] <- as.character(dtSS.table[[by]])
    dtSS.table[duplicated(interaction(dtSS.table$n,dtSS.table[[by]])), c(by) := ""]
    dtSS.table$n <- paste0("\\(n=m=\\) ",dtSS.table$n)
    dtSS.table$n[duplicated(dtSS.table$n)] <- ""
    
    setnames(dtSS.table, old = "n", new = "sample size")

    if(by=="threshold"){
        setnames(dtSS.table, old = "threshold", new = "\\(\\tau\\)")
    }
    dtSS.table[, method := factor(method, levels = c("GS","Gehan","Peron"),
                                  labels = c("Full data","Gehan", "Peron"))]
    setnames(dtSS.table, old = "method", new = "scoring rule")

    dtSS.table[, bias := gsub("<","\\(<\\)",format.pval(bias, digits = digits, eps = 10^(-digits)), fixed = TRUE)]
    dtSS.table[, empirical := format.pval(empirical, digits = digits, eps = 10^(-digits))]
    dtSS.table[, estimated := format.pval(estimated, digits = digits, eps = 10^(-digits))]
    dtSS.table[, coverage := format.pval(coverage, digits = digits, eps = 10^(-digits))]

    setnames(dtSS.table, old = "empirical", new = "empirical \\(\\sigma_{\\hat{\\Delta},\\hat{\\Delta}}\\)")
    setnames(dtSS.table, old = "estimated", new = "estimated \\(\\sigma_{\\hat{\\Delta},\\hat{\\Delta}}\\)")

    out <- xtable(dtSS.table, type = "latex", 
                  label = label,
                  caption = caption)

    addtorow <- list()
    addtorow$pos <- as.list(seq(from = n.method, to = NROW(dtSS.table), by=n.method))
    addtorow$command <- rep(c(rep("[2mm]", times = n.method-1),"[4mm]"),
                            times = length(addtorow$pos)/n.method)

    addtorow$pos <- addtorow$pos[-length(addtorow$pos)]
    addtorow$command <- addtorow$command[-length(addtorow$command)]
    mytable <- capture.output(print(out, add.to.row = addtorow, include.rownames=FALSE, include.colnames = TRUE,
                                    sanitize.colnames.function = identity,
                                    sanitize.text.function = identity,
                                    table.placement = "!h"))
    if(print){print(mytable) }

    ## ** export 
    return(invisible(list(table = mytable,
                          data = dtS.table)))
}

######################################################################
### FCT-gg.R ends here
