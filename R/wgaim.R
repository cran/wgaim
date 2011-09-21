
wgaim <- function (baseModel, ...)
UseMethod("wgaim")

wgaim.default <- function(baseModel, ...)
  stop("Currently the only supported method is \"asreml\"")

wgaim.asreml <- function (baseModel, phenoData, intervalObj, merge.by = NULL,
                          gen.type = "interval", method = "random", TypeI = 0.05,
                          attempts = 5, trace = TRUE, verboseLev = 0, ...)
{
    if (missing(phenoData))
        stop("phenoData is a required argument.")
    if (missing(intervalObj))
        stop("intervalObj is a required argument.")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (is.null(merge.by))
        stop("Need name of matching column to merge datasets.")
    if (is.null(other <- intervalObj$pheno[, merge.by]))
        stop("Genotypic data does not contain column \"", merge.by,
             "\".")
    if (is.null(phenoData[, merge.by]))
        stop("Phenotypic data does not contain column \"", merge.by,
            "\".")
    mby <- pmatch(as.character(other), as.character(phenoData[, merge.by]))
    if (all(is.na(mby)))
        stop("Names in Genotypic \"", merge.by, "\" column do not match any names in Phenotypic \"",
             merge.by, "\" column.")
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
    }
    if(gen.type == "interval")
        gdat <- lapply(intervalObj$geno, function(el) el$intval)
    else
        gdat <- lapply(intervalObj$geno, function(el) el$imputed.data)
    nint <- lapply(gdat, function(el) 1:ncol(el))
    lint <- unlist(lapply(nint, length))
    gnams <- paste("Chr", rep(names(intervalObj$geno), times = lint), unlist(nint), sep =".")
    geneticData <- cbind.data.frame(other, do.call("cbind", gdat))
    names(geneticData) <- c(merge.by, gnams)
    dnams <- names(phenoData)
    basedata <- eval(baseModel$call$data)
    bnams <- names(basedata)
    if (any(is.na(pmatch(bnams, dnams))))
        stop("Some baseModel data names do not match phenoDat names")
    whn <- unlist(lapply(basedata, is.numeric))
    whn <- whn[whn][1]
    diff <- unique(abs(basedata[[names(whn)]] - phenoData[[names(whn)]]))
    if (length(diff[!is.na(diff)]) > 1)
        stop("Phenotypic data is not in same order as baseModel data.\n Try reordering phenoData apppropriately\n")
    int.cnt <- 2:dim(geneticData)[2]

    # Is p>q?

    mD <- mergeData(phenoData, geneticData, merge.by)
    cnt <- mD$cnt
    asdata <- mD$asdata
    state <- rep(1, length(int.cnt))
    names(state) <- names(geneticData[, int.cnt])
    baseModel$call$data <- quote(asdata)
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
    }
    add.qtl <- baseModel
    class(add.qtl) <- "wgaim"
    add.qtl$call$group$ints <- cnt
    add.form <-  as.formula(paste("~ idv(grp('ints')) + ."))
    cat("\nRandom Effects QTL Model Iteration (1):\n")
    cat("=========================================\n")
    add.qtl$call$data <- quote(asdata)
    add.qtl <- updateWgaim(add.qtl, attempts = attempts, random. = add.form, ...)
    update <- FALSE
    which.i <- 1
    dmat <- data.frame(L0 = 0, L1 = 0, Statistic = 0, Pvalue = 0)
    vl <- cl <- oint <- ochr <- list()
    qtl.form <- as.formula(paste("~ . + idv(grp('qtls'))"))
    qtl <- c()
    repeat {
        if (update) {
            baseModel$call$data <- add.qtl$call$data <- quote(asdata)
            add.qtl$call$group$ints <- cnt
            if(method == "random"){
                add.qtl$call$group$qtls <- baseModel$call$group$qtls <- qtls.cnt
                cat("\nRandom Effects QTL Model Iteration (", which.i, "):\n")
                cat("========================================\n")
                baseModel <- v.modify(baseModel)
                if(which.i == 2)
                    baseModel <- updateWgaim(baseModel, attempts = attempts, random. = qtl.form, ...)
                else baseModel <- updateWgaim(baseModel, attempts = attempts, ...)
                sigma2 <- baseModel$sigma2
                if (names(baseModel$gammas.con[length(baseModel$gammas.con)]) == "Fixed") sigma2 <- 1
                cat("\nRandom Effects QTL plus Interval/Marker Model Iteration (", which.i,"):\n")
                cat("=============================================================\n")
                add.qtl <- v.modify(add.qtl)
                if(which.i == 2)
                    add.qtl <- updateWgaim	(add.qtl, attempts = attempts, random. = qtl.form, ...)
                else add.qtl <- updateWgaim	(add.qtl, attempts = attempts, ...)
                list.coefs <- add.qtl$coefficients$random
                zind <- grep("X\\.", names(list.coefs))
                list.coefs <- list.coefs[zind]
                names(list.coefs) <- gsub("qtls_","",names(list.coefs))
                cl[[which.i - 1]] <- list.coefs
                vl[[which.i - 1]] <- add.qtl$vcoeff$random[zind]
            }
            else {
                fix.form <- as.formula(paste(". ~ . +", qtl.x, sep = ""))
                cat("\nFixed Effects QTL Model Iteration (", which.i, "):\n")
                cat("========================================\n")
                baseModel <- updateWgaim(baseModel, attempts = attempts, fixed. = fix.form, ...)
                sigma2 <- baseModel$sigma2
                if (names(baseModel$gammas.con[length(baseModel$gammas.con)]) == "Fixed") sigma2 <- 1
                cat("\nRandom Effects QTL Model Iteration (", which.i, "):\n")
                cat("=========================================\n")
                add.qtl <- updateWgaim(add.qtl, attempts = attempts, fixed. = fix.form, ...)
                list.coefs <- add.qtl$coefficients$fixed
                zind <- grep("X\\.", names(list.coefs))
                cl[[which.i - 1]] <- rev(list.coefs[zind])
                vl[[which.i - 1]] <- rev(add.qtl$vcoeff$fixed[zind])
            }
        }
        baseLogL <- baseModel$loglik
        stat <- 2 * (add.qtl$loglik - baseLogL)
        dmat[which.i, ] <- c(baseLogL, add.qtl$loglik, stat, (1 - pchisq(stat, 1))/2)
        if (stat < qchisq(1 - 2 * TypeI, 1))
            break
        if (mD$p > mD$q)
            add.qtl$xtra <- mD$xtra
        add.qtl$call$data <- quote(asdata)
        if(method == "random") pick <- qtl.pick.new(add.qtl, intervalObj, gen.type, state, verboseLev)
        else pick <- qtl.pick.old(add.qtl, intervalObj, gen.type, state, verboseLev)
        state <- pick$state
        cnt <- pick$cnt
        qtl[which.i] <- pick$qtl
        oint[[which.i]] <- pick$oint
        ochr[[which.i]] <- pick$ochr
        qtl.x <- gsub("Chr\\.", "X.", qtl[which.i])
        if (is.null(add.qtl$xtra))
            phenoData[qtl.x] <- mD$asdata[qtl[which.i]] * 100
        else {
            tmp.1 <- cbind.data.frame(geneticData[, merge.by], geneticData[, qtl[which.i]])
            names(tmp.1) <- c(merge.by, qtl.x)
            phenoData <- cbind.data.frame(ord = 1:nrow(phenoData), phenoData)
            phenoData <- merge(phenoData, tmp.1, by = merge.by, all.x = TRUE, all.y = FALSE)
            phenoData <- phenoData[order(phenoData$ord), ]
            phenoData <- phenoData[, -2]
        }
        gbind <- (2:dim(geneticData)[2])[!as.logical(state)]
        gD <- geneticData[, -gbind]
        mD <- mergeData(phenoData, gD, merge.by)
        cnt <- mD$cnt
        qtls.cnt <- mD$qtls.cnt
        asdata <- mD$asdata
        which.i <- which.i + 1
        update <- TRUE
    }
    baseModel <- add.qtl
    class(baseModel) <- c("wgaim", "asreml")
    if (length(qtl)) {
        baseModel$QTL$method <- method
        baseModel$QTL$type <- gen.type
        baseModel$QTL$qtl <- qtl
        baseModel$QTL$effects <- cl[[which.i - 1]]
        baseModel$QTL$veffects <- vl[[which.i - 1]]
        baseModel$QTL$diag$dmat <- dmat
        baseModel$QTL$diag$vl <- vl
        baseModel$QTL$diag$cl <- cl
        baseModel$QTL$diag$oint <- oint
        if(method == "fixed") baseModel$QTL$diag$ochr <- ochr
    }
    data.name <- paste(as.character(baseModel$call$fixed[2]), "data", sep = ".")
    baseModel$call$data <- as.name(data.name)
    assign(data.name, asdata, envir = .GlobalEnv)
    baseModel
}

v.modify <- function(model) {
    which.term <- grep('qtls', names(model$G.param))
    if(length(which.term) > 0) {
        model$G.param <- v.init('qtls', model$G.param)
    }
    which.term <- grep('ints', names(model$G.param))
    if(length(which.term) > 0) {
        model$G.param <- v.init('ints', model$G.param)
    }
    model
}

v.init <- function(char, G.param) {
    which.term <- grep(char, names(G.param))
    var.terms <- grep('var', names(G.param[[which.term]][[1]]$initial))
    con.terms <- G.param[[which.term]][[1]]$con[var.terms] == 'B'
    if(any(con.terms)) {
        G.param[[which.term]][[1]]$con[con.terms] <- 'P'
        G.param[[which.term]][[1]]$initial[con.terms] <- 0.1
    }
    G.param
}

mergeData <- function(phenoData, geneticData, by) {

    int.cnt <- 2:dim(geneticData)[2]
    p <- length(int.cnt)
    whg <- !duplicated(phenoData[,by])
    whg <- geneticData[, by] %in% phenoData[whg,by]
    ids <- as.character(geneticData[whg, by])
    q <- length(ids)
    phenoData <- cbind(ord = 1:nrow(phenoData), phenoData)
    if(p > q) {
        mats <- geneticData[whg,int.cnt]
        tmat <- t(mats)
        xsvd <- svd(crossprod(tmat))
        xsvd.half <- xsvd.half <- t(xsvd$v %*% (t(xsvd$u) * sqrt(xsvd$d)))
        xsvd.inv <- solve(xsvd.half)
        xtra <- tmat %*% xsvd.inv
        xsvd.df <- as.data.frame(xsvd.half)
        names(xsvd.df) <- paste("Tint.", 1:q, sep = "")
        xsvd.df[[by]] <- ids
        asdata <- merge(phenoData, xsvd.df, all.x = TRUE, by = by)
        asdata <- asdata[order(asdata$ord), ]
        asdata <- asdata[, -2]
        cnt <- grep("Tint\\.", names(asdata))
        asdata[, cnt] <- asdata[, cnt]/100
    }
    else {
        xtra <- NULL
        asdata <- merge(phenoData, geneticData, by.x = by, by.y = by, all.x = TRUE, all.y = FALSE)
        asdata <- asdata[order(asdata$ord), ]
        asdata <- asdata[, -2]
        dnams <- names(asdata)
        cnt <- grep("Chr\\.", dnams)
        asdata[, cnt] <- asdata[, cnt]/100
    }
    qtls.cnt <- grep("X\\.", names(asdata))
    invisible(list(asdata=asdata, cnt=cnt, qtls.cnt = qtls.cnt, p=p, q=q, xtra=xtra))
}

updateWgaim <- function(object, ...){
  attempts <- list(...)$attempts
  class(object) <- "asreml"
  object <- update(object, ...)
  fwarn <- function(object) {
        mon <- object$monitor
        prgam <- mon[4:nrow(mon), (ncol(mon) - 2)]
        prgam[abs(prgam) < .Machine$double.eps] <- NA
        pc <- abs((mon[4:nrow(mon), (ncol(mon) - 1)] - prgam)/prgam)
        ifelse(pc[!is.na(pc)] > 0.01, TRUE, FALSE)
      }
  att <- 1
  while(!object$converge){
    object <- update(object, ...)
    att <- att + 1
    if (att > attempts)
      error.code("unstable")
  }
  att <- 1
  while (any(fwarn(object))) {
    object <- update(object, ...)
    att <- att + 1
    if (att > attempts) {
      error.code("converge")
      break
    }
  }
  object
 }

error.code <- function(ec = NULL){
    if(is.null(ec))
        stop("wgaim has been halted!")
    switch(ec,
           unstable = stop("Likelihood not converging or one or more parameters are unstable,try increasing the number of attempts or change the model"),
           converge = message("Warning message:Parameter(s) not converging: one or more parameters are unstable. Continuing with QTL analysis ...."))
}

getQTL <- function (object, intervalObj)
{
  spe <- strsplit(substring(unlist(names(object$QTL$effects)), 3),"\\.")
  wchr <- unlist(lapply(spe, function(el) el[1]))
  wint <- as.numeric(unlist(lapply(spe, function(el) el[2])))
  qtlm <- matrix(ncol = 6, nrow = length(wchr))
  for (i in 1:length(wchr)) {
    lhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
    qtlm[i, 1:4] <-  c(wchr[i], wint[i], names(lhmark), round(lhmark, 2))
    if(object$QTL$type == "interval"){
        if(length(intervalObj$geno[[wchr[i]]]$map) > 1)
            rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i] + 1]
        else rhmark <- intervalObj$geno[[wchr[i]]]$map[wint[i]]
        qtlm[i, 5:6] <- c(names(rhmark), round(rhmark, 2))
   }
   else qtlm <- qtlm[,-c(5:6)]
 }
  qtlm
}

summary.wgaim <- function (object, intervalObj, LOD = TRUE, ...) {
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (is.null(qtls <- object$QTL$effects)) {
        cat("There are no significant putative QTL's\n")
        return()
    }
    else {
        sigma2 <- object$sigma2
        if(names(object$gammas.con[length(object$gammas.con)]) == "Fixed")
            sigma2 <- 1
        zrat <- qtls/sqrt(object$QTL$veffects * sigma2)
        if(object$QTL$method == "random"){
            Pvalue <- pnorm(-zrat)
            Pvalue[qtls > 0] <- 1-Pvalue[qtls > 0]
            adds <- round(Pvalue, 3)
            addn <- "Prob"
        }
        else {
            adds <- cbind(round(zrat, 3), round(2*(1-pnorm(abs(zrat))), 3))
            addn <- c("z.ratio", "Pr(z)")
        }
        qtlmat <- as.data.frame(matrix(getQTL(object, intervalObj), nrow = length(qtls)))
        qtlmat <- cbind.data.frame(qtlmat[,c(1,3:dim(qtlmat)[2])], round(qtls, 3), adds)
        if(object$QTL$type == "interval")
            collab <- c("Chromosome", "Left Marker", "dist(cM)",
                        "Right Marker", "dist(cM)", "Size", addn)
        else
            collab <- c("Chromosome", "Marker", "dist(cM)", "Size", addn)
        if (LOD) {
            qtlmat <- cbind.data.frame(qtlmat, round(0.5 * log(exp(zrat^2),
                 base = 10), 3))
             collab <- c(collab, "LOD")
        }
        qtlmat <- qtlmat[order(qtlmat[, 1], as.numeric(as.character(qtlmat[,3]))), ]
        rownames(qtlmat) <- as.character(1:length(qtls))
        names(qtlmat) <- collab
    }
    qtlmat
}

qtlTable <- function(..., intervalObj = NULL, labels = NULL, columns = "all"){
    if(is.null(intervalObj))
        stop("Argument intervalObj cannot be NULL")
    nams <- unlist(lapply(list(...), function(el) el$QTL$type))
    if(length(unique(nams)) > 1)
        stop("Models must have been analysed with the same genetic type.")
    dots <- list(...)
    qlist <- lapply(dots, function(el, intervalObj, columns){
                    summ <- summary(el, intervalObj)
                    if(is.numeric(columns))
                        summ  <- summ[,columns]
                    summ}, intervalObj, columns)
    ql <- unlist(lapply(qlist, function(el) dim(el)[1]))
    if(!is.null(labels)){
       if(length(labels) != length(dots))
           stop("Length of labels is not equal to the number of models.")
   }
    else labels <- unlist(lapply(dots, function(el) deparse(el$call$fixed[[2]])))
    wl <- cumsum(c(1, ql[1:(length(ql) - 1)]))
    qt <- do.call("rbind.data.frame", qlist)
    qt <- cbind.data.frame(Trait = NA, qt)
    qt$Trait[wl] <- labels
    qt
}

print.wgaim <- function (x, intervalObj, ...)
{
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "interval"))
        stop("intervalObj is not of class \"interval\"")
    if (is.null(x$QTL$effects))
        cat("There are no significant putative QTL's\n")
    else {
        qtlm <- getQTL(x, intervalObj)
        for (z in 1:nrow(qtlm)) {
            int <- paste(qtlm[z, 1], qtlm[z, 2], sep = ".")
            if(x$QTL$type == "interval")
                cat("\nPutative QTL found on the interval", int,
                    "\nLeft-hand marker is", qtlm[z, 3], "\nRight-hand marker is",
                    qtlm[z, 5], "\n")
            else cat("\nPutative QTL found close to marker", int,
                    "\nMarker is", qtlm[z, 3], "\n")
        }
    }
}

tr.wgaim <- function (object, iter = 1:length(object$QTL$effects), diag.out = TRUE,
    ...)
{
    dots <- list(...)
    if (!is.na(pmatch("digits", names(dots))))
        dig <- dots$digits
    else dig <- options()$digits
    cl <- object$QTL$diag$cl
    vl <- object$QTL$diag$vl
    sigma2 <- object$sigma2
    if(names(object$gammas.con[length(object$gammas.con)]) == "Fixed")
        sigma2 <- 1
    zrl <- lapply(1:length(object$QTL$effects), function(i, cl, vl, sigma2)
                  cl[[i]]/(sqrt(vl[[i]]*sigma2)), cl = cl, vl = vl, sigma2 = sigma2)
    if (any(ret <- is.na(pmatch(iter, 1:length(zrl))))) {
        warning("\"iter\" values outside expected range .. using ones that are in iteration range")
        iter <- iter[!ret]
    }
    if(object$QTL$method == "random"){
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- pnorm(-el)
            pv[el > 0] <- 1 - pv[el > 0]
            pv <- round(pv, dig)
            c(pv, rep(NA, len - length(pv)))
    }, len = length(zrl), dig = dig)
    } else {
        pvals <- lapply(zrl, function(el, len, dig) {
            pv <- round(2 * (1 - pnorm(abs(el))), dig)
            c(pv, rep(NA, len - length(pv)))
        }, len = length(zrl), dig = dig)
    }
    qtlmat <- do.call("rbind", pvals)
    qnams <- gsub("X\\.", "", names(object$QTL$effects))
    dimnames(qtlmat) <- list(paste("Iter", 1:length(zrl), sep = "."),
        qnams)
    cat("\nIncremental QTL P-value Matrix.\n")
    cat("===============================\n")
    qtlmat[qtlmat < 0.001] <- "<0.001"
    qtlmat[is.na(qtlmat)] <- ""

    print.default(qtlmat[iter, 1:iter[length(iter)]], quote = FALSE,
        right = TRUE, ...)
    if (diag.out) {
        cat("\nOutlier Detection Diagnostic.\n")
        cat("=============================\n")
        dmat <- round(as.matrix(object$QTL$diag$dmat), dig)
        dimnames(dmat)[[1]] <- paste("Iter", 1:(length(zrl) +
            1), sep = ".")
        dmat[, 4][dmat[, 4] < 0.001] <- "<0.001"
        print.default(dmat, quote = FALSE, right = TRUE, ...)
    }
}

qtl.pick.new <- function(asr, intervalObj, gen.type, state, verboseLev) {

    ## is it (p > q) or (q > p)?
    sigma2 <- asr$sigma2
    if(names(asr$gammas.con[length(asr$gammas.con)]) == "Fixed")
        sigma2 <- 1
    if(!is.null(xtra <- asr$xtra)) {
        whq <- unlist(eval(asr$call$group$ints))
        n.int <- length(whq)
        cat(" Predict step for outlier statistics \n")
        cat("=====================================\n")

        ## use predict for random effects a and their pev

        qlev <- diag(n.int)
        qlist <- list(as.vector(qlev))
        names(qlist) <- "ints"

        pv <- predict(asr, classify = "ints",
                      only = "ints",
                      levels=qlist, vcov=TRUE, maxiter=1)

        ## covariance matrix

        gamma <- pv$gammas[grep("ints", names(asr$gammas))]

        ## atilde and var(atilde)

        atilde <- pv$predictions$pvals[, 'predicted.value']
        atilde <- xtra %*% atilde

        pev <- pv$predictions$vcov
        vatilde <- sigma2 * gamma * diag(n.int) - pev
        vatilde <- xtra %*% vatilde %*% t(xtra)
        vatilde <- diag(vatilde)
        gnams <- names(state)[as.logical(state)]
    }
    else {
        ## atilde and var(atilde) for p < q

        gamma <- asr$gammas[grep("ints", names(asr$gammas))]
        rnams <- names(asr$coefficients$random)
        grp <- grep("ints", rnams)
        atilde <- asr$coefficients$random[grp]
        pevar <- sigma2 * asr$vcoeff$random[grp]
        vatilde <- sigma2 * gamma - pevar
        gnams <- substring(rnams[grp], first=6)
    }
    names(atilde) <- gnams
    names(vatilde) <- gnams
    ntj2 <- ifelse(!is.na(atilde^2/vatilde), atilde^2/vatilde, 0)
    names(ntj2) <- gnams
    qtl <- nint <- names(ntj2)[ntj2 == max(ntj2)]
    qsp <- unlist(strsplit(nint, split="\\."))
    wint <- as.numeric(qsp[3]); wchr <- qsp[2]
    whc <- unlist(lapply(strsplit(gnams, split="\\."), function(el) el[2]))
    chri <- ntj2[whc %in% wchr]
    aint <- as.numeric(unlist(lapply(strsplit(names(chri), split="\\."), function(el) el[3])))
    if(verboseLev > 0) {
        cgen <- "Interval"
        if(gen.type == "marker") cgen <- "Marker"
        cat(cgen, "outlier statistics \n")
        cat("=============================================== \n")
        for(i in 1:length(ntj2))
            cat(cgen, names(ntj2)[i], "Outlier Statistic ", ntj2[i],"\n")
        cat("=============================================== \n\n")
   }

    ## fill out interval stats and update state

    oint <- c(ntj2)
    tint <- state
    tint[as.logical(state)] <- oint
    oint <- tint
    state[nint] <- 0
    message("Found QTL on chromosome ", wchr, " ", gen.type, " ", wint)
    list(state = state, qtl = qtl, oint = oint)
}

qtl.pick.old <- function(asr, intervalObj, gen.type, state, verboseLev) {

    ## is it (p > q) or (q > p)?
    sigma2 <- asr$sigma2
    if(names(asr$gammas.con[length(asr$gammas.con)]) == "Fixed")
        sigma2 <- 1
    if(!is.null(xtra <- asr$xtra)) {
        whq <- unlist(eval(asr$call$group))
        n.int <- length(whq)
        cat(" Predict step for outlier statistics \n")
        cat("=====================================\n")

        ## use predict for random effects a and their pev

        qlev <- diag(n.int)
        qlist <- list(as.vector(qlev))
        names(qlist) <- "ints"

        pv <- predict(asr, classify = "ints",
                      only = "ints",
                      levels=qlist, vcov=TRUE, maxiter=1)

        ## covariance matrix

        gamma <- pv$gammas[grep("ints", names(asr$gammas))]

        ## atilde and var(atilde)

        atilde <- pv$predictions$pvals[, 'predicted.value']
        atilde <- xtra %*% atilde

        pev <- pv$predictions$vcov
        vatilde <- sigma2 * gamma * diag(n.int) - pev
        vatilde <- xtra %*% vatilde %*% t(xtra)
        vatilde <- diag(vatilde)
        gnams <- names(state)[as.logical(state)]
    }
    else {
        ## atilde and var(atilde) for p < q

        gamma <- asr$gammas[grep("ints", names(asr$gammas))]
        rnams <- names(asr$coefficients$random)
        grp <- grep("ints", rnams)
        atilde <- asr$coefficients$random[grp]
        pevar <- sigma2 * asr$vcoeff$random[grp]
        vatilde <- sigma2 * gamma - pevar
        gnams <- substring(rnams[grp], first=6)
    }
    names(atilde) <- gnams
    names(vatilde) <- gnams
    ntj2 <- ifelse(!is.na(atilde^2/vatilde), atilde^2/vatilde, 0)
    names(ntj2) <- gnams

####### Chromosome first

    chr.names <- names(intervalObj$geno)
    nochr <- length(chr.names)
    nsumtj2 <- c()
    for(c in 1:nochr){
        wha <- grep(chr.names[c], gnams)
        whc <- unlist(lapply(strsplit((gnams)[wha], split='\\.'), function(x) x[2]))
        whl <- whc %in% chr.names[c]
        wha <- wha[whl]
        catilde <- atilde[wha]
        nums <- catilde * catilde
        names(nums) <- names(atilde)[wha]
        dens <- vatilde[wha]
        nsumtj2[c] <- ifelse(!is.na(sum(nums)/sum(dens)),sum(nums)/sum(dens),0)
    }
    names(nsumtj2) <- chr.names # names(nmark)
    ochr <- c(nsumtj2)
    chind <- (1:nochr)[nsumtj2 == max(nsumtj2)]
    wchr <- chr.names[chind] # names(nmark)[chind]
    chri <- ntj2[grep(wchr, gnams)]
    whc <- unlist(lapply(strsplit(names(chri), split='\\.'), function(x) x[2]))
    whl <- whc %in% wchr
    chri <- chri[whl]
    wint <- (1:length(chri))[chri == max(chri)]
    qtl <- nint <- names(chri)[wint]
   if(verboseLev > 0) {
     cat("\n Selection of chromosome using the AOM statistic\n")
     cat("=============================================== \n")
     for(i in 1:nochr)
       cat(" Chromosome ", chr.names[i], "Outlier Statistic ", nsumtj2[i], "\n")
     cat("============================================= \n\n")
     cgen <- "Interval"
     if(gen.type == "marker") cgen <- "Marker"
     cat(cgen, "outlier statistics \n")
     cat("=============================================== \n")
     for(i in 1:length(chri))
         cat(cgen, names(chri)[i], "Outlier Statistic ", chri[i],"\n")
     cat("=============================================== \n\n")
   }

    ## fill out interval stats and update state
    oint <- c(ntj2)
    tint <- state
    tint[as.logical(state)] <- oint
    oint <- tint
    state[nint] <- 0
    message("Found QTL on chromosome ", wchr, " ", gen.type, " ", wint)
    list(state = state, qtl = qtl, ochr = ochr, oint = oint)
}

cross2int <- function(fullgeno, missgeno = "MartinezCurnow", rem.mark = TRUE, id = "id",
                      subset = NULL){
  if(!inherits(fullgeno, "bc"))
    stop("This function is restricted to populations containing only two genotypes.")
  fullgeno <- drop.nullmarkers(fullgeno)
  if(!(id %in% names(fullgeno$pheno)))
    stop("The unique identifier for the genotypic rows, ", deparse(substitute(id)), ",cannot be found in genotypic data")
  if(!is.null(subset))
    fullgeno <- subset(fullgeno, chr = subset)
  if(rem.mark) {
    tpheno <- fullgeno$pheno
    fullgeno <- fix.map(fullgeno, getwd(), "dummy", id)
    fullgeno$pheno <- tpheno
  }
  lid <- as.character(fullgeno$pheno[[id]])
  mtype <- c("Broman", "MartinezCurnow")
  if(is.na(type <- pmatch(missgeno, mtype)))
      stop("Missing marker type must be one of \"Broman\" or \"MartinezCurnow\". Partial matching is allowed.")
  missgeno <- mtype[type]
  if(missgeno == "Broman")
      fullgeno <- argmax.geno(fullgeno)
  fullgeno$geno <- lapply(fullgeno$geno, function(el, missgeno, lid){
      row.names(el$data) <- as.character(lid)
      if(length(el$map) == 1){
          el$dist <- 0; el$theta <- 0; elambda <- 1/2
          names(el$dist) <- names(el$map)
          el$data[is.na(el$data)] <- 0
          el$imputed.data <- el$data
          el$imputed.data[el$imputed.data == 2] <- - 1
          el$intval <- as.matrix(el$imputed.data/2, ncol = 1)
          dimnames(el$intval)[[2]] <- names(el$map)
      }
      else {
          el$dist <- diff(el$map)/100
          el$theta <- 0.5*(1-exp(-2*el$dist))
          elambda <- el$theta/(2*el$dist*(1-el$theta))
          if(missgeno == "MartinezCurnow"){
              el$imputed.data <- el$data
              el$imputed.data[el$imputed.data == 2] <- -1
              el$imputed.data <- miss.q(el$theta, el$imputed.data)
          }
          else {
              el$imputed.data <- el$argmax
              el$imputed.data[el$imputed.data == 2] <- -1
              dimnames(el$imputed.data)[[1]] <- dimnames(el$data)[[1]]
          }
          lambda <- addiag(elambda,-1) + addiag(c(elambda,0),0)
          lambda <- lambda[,-dim(lambda)[2]]
          el$intval <- el$imputed.data %*% lambda
          dimnames(el$intval)[[2]] <- names(el$dist)
      }
      el
  }, missgeno, lid)
  class(fullgeno) <- c(class(fullgeno), "interval")
  fullgeno
}

fix.map <- function(full.data, dir, filename, id){
    full.data <- est.rf(full.data)
    mymat <- full.data$rf
    mymat[!lower.tri(mymat)] <- NA
    mylist <- which(mymat==0,arr.ind = TRUE)
    if(length(mylist)==0)
        return(full.data)
    else
      mylist <- mylist[,2]
    index <- 1:sum(nmar(full.data))
    indexnames <- rep(names(cumsum(nmar(full.data))),nmar(full.data))
    mymat <- as.data.frame(mymat)
    finallist <- cbind(names(mylist),names(mymat)[mylist],indexnames[mylist])
    newmap <- drop.markers(full.data, unique(finallist[,1]))
    cat("\nDropping coincident markers.....\n")
    newdat <- do.call("rbind", lapply(newmap$geno, function(el) t(el$data)))
    newdat <- cbind.data.frame(row.names(newdat), rep(names(nmar(newmap)), times = nmar(newmap)), newdat)
    names(newdat) <- c(id, "", as.character(newmap$pheno[,id]))
    fsub <- paste(filename, ".csv", sep ="")
    cat("Creating new map",deparse(substitute(fsub)),"....\n\n")
    write.csv(newdat, paste(dir, fsub, sep="\\"), row.names=FALSE)
    # Read newmap back in.
    full.data <- read.cross("csvr", dir, file=fsub, genotypes=c("1","2"), na.strings="NA")
    full.data$cor.markers <- finallist
    full.data
  }

addiag <- function(x = 1, di = 0, nrow.arg, ncol.arg = n)
{
    if(is.matrix(x)) {
        k <- ifelse(col(x) == (row(x) + di), TRUE, FALSE)
        return(x[k])
    }
    if(missing(x))
        n <- nrow.arg
    else if(length(x) == 1 && missing(di) && missing(nrow.arg) && missing(ncol.arg)) {
        n <- as.integer(x)
        x <- 1
    }
    else n <- length(x)
    if(!missing(nrow.arg))
        n <- nrow.arg
    k <- abs(di)
    p <- ncol.arg + k
    n <- n + k
    m <- matrix(0, n, p)
    k <- ifelse(col(m) == (row(m) + di), TRUE, FALSE)
    m[k] <- x
    m
}

miss.q <- function(theta, chr){
    n <- ncol(chr)
    m <- nrow(chr)
    mk <- (1:m)[apply(chr, 1, function(el) any(is.na(el)))]
    for(k in mk){
      whk <- (1:n)[is.na(chr[k,])]
      if(length(whk) == n) {
          chr[k, ] <- 0
          warning("Line ", dimnames(chr)[[1]][k], " has missing values across the whole of a chromosome, .. These have been replaced by 0's.")
      }
      else {
          if((n %in% whk)) {
              if(length(whk) >= 2)
                  whk <- c(whk[length(whk)], whk[-length(whk)])
          }
          for(i in whk){
              thetaR <- thetaL <- 0
              if(i == n){
                  for(l in 1:(n - 1)) {
                      thetaL <- thetaL + theta[n - l] - 2*thetaL*theta[n - l]
                      xL <- chr[k, n - l]
                      if(!is.na(xL))
                          break
                  }
                  chr[k, i] <- as.numeric(xL)*(1 - 2*thetaL)
              }
              else {
                  for(j in (i + 1):n) {
                      thetaR <- thetaR + theta[j - 1] - 2*thetaR*theta[j - 1]
                      xR <- chr[k, j]
                      if(!is.na(xR))
                          break
                  }
                  if(i == 1)
                      chr[k, i] <- as.numeric(xR)*(1 - 2*thetaR)
                  else {
                      thetaL <- theta[i - 1]
                      xL <- chr[k, i - 1]
                      thetaLR <- thetaL + thetaR - 2*thetaL*thetaR
                      if(thetaLR == 0)
                          chr[k, i] <- xL
                      else {
                          lambda <- (thetaR*(1 - thetaR)*(1 - 2*thetaL))/(thetaLR*(1 - thetaLR))
                          rho <- (thetaL*(1 - thetaL)*(1 - 2*thetaR))/(thetaLR*(1 - thetaLR))
                          chr[k, i] <- as.numeric(xL)*lambda + as.numeric(xR)*rho
            }
          }
        }
      }
    }
  }
    chr
  }

link.map <- function(object, ...)
  UseMethod("link.map")

link.map.cross <- function(object, chr, max.dist, marker.names = "markers", tick = FALSE, squash = TRUE, m.cex = 0.6, ...){
  circ <- function(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
    ((x - shiftx)^2)/elx + ((y - shifty)^2)/ely
  dots <- list(...)
  old.xpd <- par("xpd")
  par(xpd = TRUE)
  on.exit(par(xpd = old.xpd))
  map <- pull.map(object)
  if(!missing(chr)) {
    if(any(is.na(pmatch(chr, names(map)))))
      stop("Some names of chromosome(s) subset do not match names of map.")
    map <- map[chr]
  }
  n.chr <- length(map)
  mt <- list()
  if(!missing(max.dist))
    map <- lapply(map, function(el, max.dist) el[el < max.dist], max.dist)
  maxlen <- max(unlist(lapply(map, max)))
  if(is.null(marker.names)) {
    chrpos <- 1:n.chr
    thelim <- range(chrpos) + c(-0.5, 0.5)
  }
  else {
    if(all(is.na(pmatch(marker.names, c("markers","dist")))))
      stop("marker.names argument must be either \"dist\" or \"markers\"")
    if(!is.na(pmatch("cex", names(dots))))
      dots$cex <- NULL
#    else cex <- par("cex")
    if(!squash)
      chrpos <- seq(1, n.chr * 3, by = 3)
    else
      chrpos <- seq(1, n.chr * 2, by = 2)
    thelim <- range(chrpos) + c(-1.6, 1.35)
    for(i in 1:n.chr) {
        mt[[i]] <- map[[i]]
      if(length(mt[[i]]) > 1){
          conv <- par("pin")[2]/maxlen
      for(j in 1:(length(mt[[i]]) - 1)){
          ch <- mt[[i]][j + 1]*conv - (mt[[i]][j]*conv + 10*par("csi")*m.cex/9)
        if(ch < 0){
            temp <- mt[[i]][j + 1]*conv + abs(ch)
          mt[[i]][j + 1] <- temp/conv
        }
      }
      }
    }
    maxlen <- max(unlist(lapply(mt, max)))
    names(mt) <- names(map)
  }
  plot(0, 0, type = "n", ylim = c(maxlen, 0), xlim = thelim,
       xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome",
       axes = FALSE, ...)
  axis(side = 2,  ylim = c(maxlen, 0))
  pins <- par()$plt
  for(i in 1:n.chr) {
    if(!is.null(marker.names)) {
       if(marker.names == "dist")
         alis <- list(x = chrpos[i] + 0.50, y =  mt[[i]], labels = as.character(round(map[[i]], 2)),
                      adj = c(0, 0.5), cex = m.cex)
       else
         alis <- list(x = chrpos[i] + 0.50, y = mt[[i]], labels = names(map[[i]]), adj = c(0, 0.5), cex = m.cex)
      do.call("text", c(alis, dots))
      segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 0.3, map[[i]])
      segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 0.4, mt[[i]])
      segments(chrpos[i] + 0.40, mt[[i]], chrpos[i] + 0.45, mt[[i]])
    }
    barl <- chrpos[i]- 0.03
    barr <- chrpos[i]+ 0.03
    segments(barl, min(map[[i]]), barl, max(map[[i]]), lwd = 1)
    segments(barr, min(map[[i]]), barr, max(map[[i]]), lwd = 1)
    segments(barl - 0.17, map[[i]], barr + 0.17, map[[i]])
    # attempt to put curves at ends of chromosomes
    xseq <- seq(barl, barr, length = 20) - chrpos[i]
    yseq <- circ(xseq, xseq, ely = 1, elx = 0.07/maxlen)
    yseq <- yseq - max(yseq)
    lines(xseq + chrpos[i], min(map[[i]]) + yseq)
    lines(xseq + chrpos[i], max(map[[i]]) - yseq)
  }
  axis(side = 1, at = chrpos, labels = names(map), tick = FALSE)
  if(is.na(pmatch("main", names(dots))) & !as.logical(sys.parent()))
    title("Genetic Map")
  invisible(list(mt = mt, map = map, chrpos = chrpos))
}

#savePlot("qtlplot", "pdf")

link.map.wgaim <- function (object, intervalObj, chr, max.dist, marker.names = "markers",
    list.col = list(q.col = "light blue", m.col = "red", t.col = "light blue"),
    list.cex = list(t.cex = 0.6, m.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    dots <- list(...)
    if (missing(intervalObj))
        stop("intervalObj is a required argument")
    if (!inherits(intervalObj, "cross"))
        stop("intervalObj is not of class \"cross\"")
    if (!length(wchr <- object$QTL$effects)) {
        warning("There are no significant QTL's. Plotting map only...")
        link.map(intervalObj, chr, max.dist, marker.names = marker.names,
            tick = tick, squash = FALSE, ...)
        return(invisible())
    }
    qtlm <- getQTL(object, intervalObj)
    wchr <- qtlm[,1]
    if(object$QTL$type == "interval")
        qtlm <- qtlm[,3:6]
    else
        qtlm <- cbind(qtlm[,3:4],qtlm[,3:4])
    if (is.null(list.col$q.col))
        list.col$q.col <- "light blue"
    if (is.null(list.col$t.col))
        list.col$t.col <- list.col$q.col
    if (missing(chr))
        chr <- unique(wchr)[order(unique(wchr))]
    if(is.null(list.cex$m.cex))
        list.cex$m.cex <- 0.6
    if(is.null(list.cex$t.cex))
        list.cex$t.cex <- 0.6
    lmap <- link.map(intervalObj, chr, max.dist, marker.names = marker.names,
        tick = tick, squash = TRUE, m.cex = list.cex$m.cex, ...)
    map <- lmap$map
    if (is.null(trait <- trait.labels))
        trait <- rep(as.character(object$call$fixed[[2]]), length(wchr))
    if (length(trait.labels) == 1)
        trait <- rep(trait.labels, length(wchr))
    qtrait <- unique(trait)
    qtlm <- cbind.data.frame(qtlm, trait = factor(trait, levels = unique(trait)))
    qtlm[, 2] <- as.numeric(as.character(qtlm[, 2]))
    qtlm[, 4] <- as.numeric(as.character(qtlm[, 4]))
    qtlList <- lapply(split(qtlm, wchr), function(el) el[order(el[,
        2]), ])
    qtlm <- do.call("rbind", qtlList)
    wchr <- wchr[order(wchr)]
    if (!missing(chr)) {
        if (any(is.na(wh <- pmatch(wchr, chr, dup = TRUE)))) {
            warning("Some QTL's exist outside chromosome(s) subset, Omitting QTL's....")
            qtlm <- qtlm[!is.na(wh), ]
            wchr <- wchr[!is.na(wh)]
        }
    }
    if (!missing(max.dist)) {
        rml <- qtlm[, 4] > max.dist
        if (any(rml)) {
            warning("Some QTL regions outside maximum distance specified. Omitting QTL's....")
            qtlm <- qtlm[!rml, ]
            wchr <- wchr[!rml]
        }
    }
    n.chr <- length(map)
    maxlen <- max(unlist(lapply(map, max)))
    chrpos <- lmap$chrpos
    mt <- lmap$mt
    if (!is.na(cind <- pmatch("col", names(dots))))
        dots <- dots[-cind]
    if (is.null(dim(qtlm)))
        qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)
    tlis <- list()
    if (!is.na(pmatch("cex", names(dots))))
        p.cex <- dots$cex
    else p.cex <- par("cex")
    dots$cex <- NULL
    for (i in 1:n.chr) {
        conv <- par("pin")[2]/maxlen
        ind <- wchr %in% names(map)[i]
#        if (as.logical(length(ind <- wchr %in% names(map)[i]))) {
         if(any(ind)){
            tlis[[i]] <- as.vector(qtlm[ind, 2] + qtlm[ind, 4])/2
            names(tlis[[i]]) <- as.character(qtlm[ind, 5])
            if (length(tlis[[i]]) > 1) {
                for (j in 1:(length(tlis[[i]]) - 1)) {
                  ch <- tlis[[i]][j + 1] * conv - (tlis[[i]][j] *
                    conv + 10 * par("csi") * list.cex$t.cex/9)
                  if (ch < 0) {
                    temp <- tlis[[i]][j + 1] * conv + abs(ch)
                    tlis[[i]][j + 1] <- temp/conv
                  }
                }
            }
        }
    }
    qtld <- qtlm[, 1:4]
    nodup <- !duplicated(do.call("paste", qtld))
    qtls <- qtld[nodup, ]
    whd <- pmatch(do.call("paste", qtld), do.call("paste", qtls),
        dup = TRUE)
    dlis <- split(as.character(qtlm[, 5]), whd)
    qtlm <- qtls
    wchr <- wchr[nodup]
    for (i in 1:n.chr) {
        ind <- wchr %in% names(map)[i]
        if(any(ind)){
            ind <- (1:length(wchr))[ind]
            for (j in ind) {
                if (!is.null(marker.names)) {
                    wh <- mt[[i]][pmatch(c(as.character(qtlm[j,
                      1]), as.character(qtlm[j, 3])), names(map[[i]]))]
                  if (marker.names == "dist") {
                    dist <- map[[i]][pmatch(c(as.character(qtlm[j,
                      1]), as.character(qtlm[j, 3])), names(map[[i]]))]
                    alis <- list(x = chrpos[i] + 0.5, y = wh,
                      labels = as.character(round(as.numeric(dist),
                        2)), adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  }
                  else alis <- list(x = chrpos[i] + 0.5, y = wh,
                    labels = names(wh), adj = c(0, 0.5), col = list.col$m.col, cex = list.cex$m.cex)
                  do.call("text", c(alis, dots))
                }
                yv <- c(qtlm[j, 2], qtlm[j, 4])
                yv <- c(yv, rev(yv))
                dind <- dlis[[j]]
                q.cols <- list.col$q.col[pmatch(dind, qtrait)]
                qind <- 1:length(dind)
                if (length(dlis[[j]]) > 1) {
                  int <- seq(chrpos[i] - 0.2, chrpos[i] + 0.2,
                    length = length(dind) + 1)
                  for (k in 1:length(dind)) {
                    xv <- c(rep(int[k], 2), rep(int[k + 1], 2))
                    if(object$QTL$type == "interval")
                        polygon(xv, y = yv, border = NA, col = q.cols[k])
                    else {
                        plis <- list(x = (xv[1] + xv[3])/2,y = yv[1], col = q.cols[k], cex = p.cex)
                        do.call("points", c(plis, dots))
                    }
                  }
                }
                else {
                  xv <- c(rep(chrpos[i] - 0.2, 2), rep(chrpos[i] +
                    0.2, 2))
                  if(object$QTL$type == "interval")
                      polygon(xv, y = yv, border = NA, col = q.cols)
                  else {
                      plis <- list(x = chrpos[i], y = yv[1], col = q.cols, cex = p.cex)
                      do.call("points", c(plis,dots))
                  }
                }
                segments(chrpos[i] - 0.25, yv[1], chrpos[i] -
                  0.25, yv[2])
                segments(chrpos[i] - 0.25, sum(yv[1:2])/2, chrpos[i] -
                  0.3, sum(yv[1:2])/2)
                segments(chrpos[i] - 0.3, sum(yv[1:2])/2, chrpos[i] -
                  0.4, tlis[[i]][qind])
                segments(chrpos[i] - 0.4, tlis[[i]][qind], chrpos[i] -
                  0.45, tlis[[i]][qind])
                if (length(list.col$t.col) > 1)
                  t.cols <- list.col$t.col[pmatch(dind, qtrait)]
                else t.cols <- list.col$t.col
                text(chrpos[i] - 0.5, tlis[[i]][qind], names(tlis[[i]][qind]),
                  adj = c(1, 0.3), col = t.cols, cex = list.cex$t.cex)
                tlis[[i]] <- tlis[[i]][-qind]
            }
        }
        segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2,
            map[[i]])
    }
    if (is.na(pmatch("main", names(dots))))
        title("Genetic Map with QTL")
}

link.map.default <- function (object, intervalObj, chr, max.dist, marker.names = "markers",
    list.col = list(q.col = rainbow(length(object)), m.col = "red", t.col = rainbow(length(object))),
    list.cex = list(m.cex = 0.6, t.cex = 0.6), trait.labels = NULL, tick = FALSE, ...)
{
    old.par <- par(no.readonly = TRUE)
    par(mar = c(5, 4, 5, 2) + 0.1)
    par(xpd = TRUE)
    on.exit(par(old.par))
    dlist <- list()
    lclass <- lapply(object, function(el) {
        if (!inherits(el, "wgaim"))
            stop("list objects need to inherit from class \"wgaim\"")
        class(el)[1]
    })
    lclass <- unique(unlist(lclass))
    if (any(is.na(pmatch(lclass, "wgaim"))))
        stop("link.map method is only for list objects of class \"wgaim\"")
    type <- unique(unlist(lapply(object, function(el) el$QTL$type)))
    if(length(type) > 1)
        stop("Models need to have same analyses \"type\".")
    dlist$QTL$type <- type
    effects <- lapply(object, function(el) el$QTL$effects)
    len <- lapply(effects, length)
    if (!is.null(trait.labels)) {
        if (length(trait.labels) != length(object))
            stop("Length of trait labels does not equal number of models specified.")
        trait <- rep(trait.labels, times = len)
    }
    else {
        trait <- unlist(lapply(object, function(el) rep(as.character(el$call$fixed[[2]]),
            length(el$QTL$effects))))
    }
    dlist$QTL$effects <- unlist(effects)
    class(dlist) <- "wgaim"
    if (length(list.col$q.col) != length(object)) {
        warning("QTL colours not of the same length as the number of traits,\n Choosing \"q.col = rainbow(",
            length(object), ")\".")
        list.col$q.col <- rainbow(length(object))
    }
    if (length(list.col$t.col) != length(object)) {
        if (length(list.col$t.col) != 1) {
            warning("Inappropriate length for trait name colours, using QTL colours")
            list.col$t.col <- list.col$q.col
        }
    }
    link.map(dlist, intervalObj, chr, max.dist, marker.names = marker.names,
        list.col = list.col, list.cex = list.cex, trait.labels = trait, tick = tick, ...)
  }

out.stat <- function (object, intervalObj, int = TRUE, iter = NULL, chr = NULL,
    ...)
{
    mod <- function(x, num) {
        x <- x/num
        val <- num * (x - floor(x))
        val[val == 0] <- num
        val
    }
    dots <- list(...)
    cex <- dots$cex
    if (is.null(qtl <- object$QTL$effects))
        stop("There are no outlier statistics to plot.")
    qtl <- gsub("X.","",names(qtl))
    if(object$QTL$type == "interval")
        dist <- lapply(intervalObj$geno, function(el) {
                        if(length(el$map) == 1) {
                            el$dist <- 0.05
                            names(el$dist) <- names(el$map)
                            el$dist
                        }
                        else el$dist })
    else
        dist <- lapply(intervalObj$geno, function(el) {
                           tel <- c(0.05, el$dist)
                           names(tel)[1] <- names(el$map)[1]
                           tel })
    ln <- length(qtl)
    ochr <- object$QTL$diag$ochr
    oint <- object$QTL$diag$oint
    cnam <- names(ochr[[1]])
    cnam <- gsub("Chr.", "", cnam)
    inam <- names(oint[[1]])
    inam <- strsplit(gsub("Chr.", "", inam), "\\.")
    chn <- unlist(lapply(inam, function(el) el[1]))
    inn <- unlist(lapply(inam, function(el) el[2]))
    if (!is.null(iter)) {
        if (is.numeric(iter))
            iter <- as.integer(iter)
        if (any(iter > ln))
            stop("\"iter\" is greater than the number of QTL found in the analysis.")
        oint <- oint[iter]
        ochr <- ochr[iter]
        ln <- length(iter)
    }
    else iter <- 1:ln
    oint <- unlist(oint)
    dint <- cbind.data.frame(int = oint, chn = rep(chn, length(iter)),
        inn = rep(inn, length(iter)))
    if (!int){
        if(object$QTL$method == "random")
            stop("This method does not produce any chromosome statistics.")
        chr <- NULL
    }
    if (!is.null(chr)) {
        dint <- dint[as.character(dint$chn) %in% chr, ]
        dist <- unlist(dist[chr])
    }
    else dist <- unlist(dist)
    labs <- unique(paste(dint$chn, ".", dint$inn, sep = ""))
    if(object$QTL$type == "interval")
        dist <- cumsum(dist) - dist/2
    else dist <- cumsum(dist)
    dint$dist <- rep(dist, length(iter))
    dint$itn <- rep(1:length(iter), each = length(labs))
    its <- paste("Iteration: ", rep(iter, each = length(labs)),
        sep = "")
    dint$its <- factor(its, levels = unique(its))
    dint$step <- rep(1:length(labs), length(iter))
    qtln <- paste(qtl[iter], unique(dint$itn), sep = ":")
    td <- dint[paste(paste(dint$chn, dint$inn, sep = "."), dint$itn,
        sep = ":") %in% qtln, ]
    prow <- length(unique(dint$its))
    if (prow > 5)
        prow <- 5
    if(int){
        slabs <- strsplit(labs, "\\.")
        clabs <- unlist(lapply(slabs, function(el) el[1]))
        tabc <- table(clabs)
        cumc <- c(1, cumsum(tabc) + 1)
        labs[-cumc] <- " "
        labs <- as.character(labs)
        print(xyplot(int ~ dist | its, type = "l", data = dint,
                     groups = chn, panel = panel.superpose, panel.groups = function(x,
                                                            y, ...) {
                         panel.xyplot(x, y, ...)
                     }, ylab = "outlier statistic", scales = list(x = list(labels = labs,
                     rot = 45, at = unique(dist), cex = 0.6)), xlab = object$QTL$type,
                     layout = c(1, prow), ...))
    } else {
        dups <- !duplicated(paste(dint[, "chn"], dint[, c("its")],
                                  sep = "."))
        dchr <- dint[dups, ]
        dchr$int <- unlist(ochr)
        names(dchr)[1] <- "chr"
        bc <- print(barchart(chr ~ chn | its, data = dchr, ylab = "outlier statistic",
            xlab = "Chromosome", layout = c(1, prow), ...))
        wchr <- dchr$chr[paste(dchr$chn, dchr$itn, sep = ":") %in%
            paste(td$chn, td$itn, sep = ":")]
        td$int <- wchr
        td$whc <- pmatch(td$chn, as.character(unique(dchr$chn)),
            duplicates.ok = TRUE)
        names(td)[1] <- "chr"
    }
    panel.locs <- as.vector(trellis.currentLayout())
    panel.locs <- panel.locs[panel.locs > 0]
    rows <- round(mod(as.vector(panel.locs), 5), 1)
    k <- 1
    for (wh in panel.locs) {
        trellis.focus("panel", row = rows[k], column = 1, highlight = TRUE)
        tq <- paste(td$chn[td$itn == wh], td$inn[td$itn == wh],
            sep = ".")
        if (int)
            panel.text(td$dist[td$itn == wh], td$int[td$itn ==
                wh], label = tq, col = 1, cex = cex)
        else {
            lims <- bc$y.limits[2] - bc$y.limits[1]
            panel.text(td$whc[td$itn == wh], td$chr[td$itn ==
                wh] + 0.1 * lims, label = tq, col = 1, cex = cex)
        }
        k <- k + 1
      }
}








