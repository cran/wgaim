wgaim <- function(baseModel, ...)
  UseMethod("wgaim")

wgaim.default <- function(baseModel, ...) 
  stop("Currently the only supported method is \"asreml\"")
  
wgaim.asreml <- function(baseModel, parentData, TypeI = 0.05, attempts = 5, trace = TRUE, ...){

  ## check data 

  if(missing(parentData))
    stop("parentData is a required argument.")
  if(!inherits(parentData, "interval"))
    stop("parentData is not of class \"interval\"")
  if(is.null(parentData$full.data))
    stop("parentData should have component \"full.data\" (see ?wmerge).")
  if(is.character(trace)){
    ftrace <- file(trace, "w")
    sink(trace, type = "output", append = FALSE)   
    on.exit(sink(type = "output"))
    on.exit(close(ftrace), add = TRUE)
  }
  dlist <- match.call(expand.dots = FALSE)$...
  asdata <- parentData$full.data
  dnams <- names(asdata)
  cnt <- grep("Chr\\.", dnams)
  asdata[,cnt] <- asdata[,cnt]/100
  state <- rep(1,length(cnt))
  names(state) <- names(asdata[,cnt])
  
  ## A few checks

  basedata <- eval(baseModel$call$data)
  bnams <- names(basedata)
  if(any(is.na(pmatch(bnams, dnams))))
    stop("Some baseModel data names do not match phenotypic parentData names")
  whn <- unlist(lapply(basedata, is.numeric))
  whn <- whn[whn][1]
  diff <- unique(abs(basedata[[names(whn)]] - asdata[[names(whn)]]))
  if(length(diff[!is.na(diff)]) > 1)
    stop("Phenotypic parentData is not in same order as baseModel data.\n Try reordering the phenotypic parentData apppropriately\n")
  
  baseModel$call$data <- quote(asdata)

  ## Check base model characteristics 
  
  if(!baseModel$converge) {
    cat("Warning: Base model has not converged. Updating base model\n")
    baseModel <- update(baseModel)
  }

  ## add ... to call
  
  if (length(dlist) > 0) {
    what <- !is.na(match(names(dlist), names(baseModel$call)))
    for (z in names(dlist)[what]) baseModel$call[[z]] <- dlist[[z]]
    if (any(!(what))) {
      baseModel$call <- c(as.list(baseModel$call), dlist[!what])
      baseModel$call <- as.call(baseModel$call)
   }
  }

  add.qtl <- baseModel

  ## Set up random formula and group structure for add.qtl 

  if(is.null(add.qtl$call$random))
    add.qtl$call$random <- as.formula(paste("~", "idv(grp('qtls'))", sep=""))    
  else
    add.qtl$call$random <- as.formula(paste("~ idv(grp('qtls'))", deparse(add.qtl$call$random[[2]], width.cutoff = 500), sep=" + "))                    
  if(is.null(add.qtl$call$group))
    add.qtl$call$group <- list("qtls" = cnt)
  else 
    add.qtl$call$group <- c(add.qtl$control$group, list('qtls' = cnt))
  
 ## initialize variables

  add.qtl$call$G.param <- baseModel$G.param
  add.qtl$call$R.param <- baseModel$R.param

  update <- FALSE
  qtl.fix <- wchr <- wint <- NULL
  i <- 1;
  gamma <- 0.1
  nochr <- nchr(parentData)
  nmark <- nmar(parentData) - 1
  names(nmark) <- paste("Chr.", names(nmark), sep ="")
  dmat <- data.frame(L0 = 0, L1 = 0, Statistic = 0, Pvalue = 0)
  zrl <- zcl <- list()
  oint <- ochr <- list()
  
  
  ## initial QTL fit

  message("Searching for QTLs......")
  cat("\nRandom Effects QTL Model Iteration (",i,"):\n")
  cat("=========================================\n")
  add.qtl <- eval(add.qtl$call)
  rnams <- names(add.qtl$coefficients$random)

  ## loop fitting and find QTLs

  fwarn <- function(model){
      mon <- model$monitor
      prgam <- mon[4:nrow(mon),(ncol(mon) - 2)]
      prgam[abs(prgam) < .Machine$double.eps] <- NA
      pc <- abs((mon[4:nrow(mon),(ncol(mon) - 1)] - prgam)/prgam)
      ifelse(pc[!is.na(pc)] > 0.01, TRUE, FALSE)
    }
      
  repeat {
    if(update) {
      cat("\nFixed Effects QTL Model Iteration (",i,"):\n")
      cat("========================================\n")
      baseModel <- eval(baseModel$call)
      b <- 1
      while(!baseModel$converge) {
        baseModel <- update.asreml(baseModel)
        b <- b + 1
        if(b > attempts) {
        message("Error message:\nLikelihood not converging or one or more parameters are unstable,
try increasing the number of attempts or change the model. For diagnostic purposes the current fixed
QTL model is returned and full details of the current working random QTL model can be found in ..$QTL")
        assign("asdata", asdata, envir = .GlobalEnv) 
        baseModel$QTL$qtlModel <- add.qtl        
        return(invisible(baseModel))
      }
      }
      zrat <- round(baseModel$coefficients$fixed/sqrt(baseModel$vcoeff$fixed*baseModel$sigma2), 2)
      zind <- grep("XChr\\.", names(zrat))
      zrl[[i - 1]] <- rev(zrat[zind])
      zcl[[i - 1]] <- rev(baseModel$coefficients$fixed[zind])
      b <- 1
      while(any(fwarn(baseModel))){
        baseModel <- update.asreml(baseModel)
        b <- b + 1
        if (b > attempts) {
          message("Warning message:\nParameter(s) not converging: one or more parameters are unstable. Continuing with QTL analysis ....\n")
          assign("asdata", asdata, envir = .GlobalEnv)
          baseModel$QTL$qtlModel <- add.qtl
          break
        }
      }
      cat("\nRandom Effects QTL Model Iteration (",i,"):\n")
      cat("=========================================\n")
      add.qtl <- eval(add.qtl$call)
    }
    baseLogL <- baseModel$loglik
    q <- 1
    while(!add.qtl$converge){
      add.qtl <- update.asreml(add.qtl)  
      q <- q + 1
      if(q > attempts) {
        message("Error message:\nLikelihood not converging or one or more parameters are unstable,
try increasing the number of attempts or change the model.  For diagnostic purposes the current fixed
QTL model is returned and full details of the current working random QTL model can be found in ..$QTL")
        assign("asdata", asdata, envir = .GlobalEnv) 
        baseModel$QTL$qtlModel <- add.qtl
        return(invisible(baseModel))
      }
    }
    q <- 1
    while(any(fwarn(add.qtl))) {
       add.qtl <- update.asreml(add.qtl)
       q <- q + 1
       if (q > attempts) {
         message("Error message:\nParameter(s) not converging: one or more parameters are unstable. Continuing with QTL analysis ....\n")
         assign("asdata", asdata, envir = .GlobalEnv)
         add.qtl$QTL$qtlModel <- add.qtl
         break                               
       }
     }
    stat <- 2*(add.qtl$loglik - baseLogL)
    dmat[i,] <- c(baseLogL, add.qtl$loglik, stat, (1 - pchisq(stat, 1))/2)
    if(stat < qchisq(1-2*TypeI,1))
      break
    
    ## Find the most likely chromosome and then the most
    ## likely interval on that chromosome
    ## that contains a putative QTL.

    gamma <- add.qtl$gammas[grep("qtls", names(add.qtl$gammas))]
    sigma2 <- add.qtl$sigma2

    rnams <- names(add.qtl$coefficients$random)
    grp <- grep("qtls", rnams)
    bnams <-  rnams[grp]
    blups <- add.qtl$coefficients$random[grp]
    blups[!as.logical(state)] <- 0
    pevar <- sigma2*add.qtl$vcoeff$random[grp]
    varu <- sigma2*gamma - pevar
    oint[[i]] <- ntj2 <- blups^2/varu
    nsumtj2 <- c()
    for(c in 1:nochr){
      wh.m <- grep(names(nmark)[c], bnams)
      nsumtj2[c] <- sum(blups[wh.m]^2)/sum(varu[wh.m])
    }
    names(nsumtj2) <- names(nmark)
    ochr[[i]] <- nsumtj2
    chind <- (1:nochr)[nsumtj2 == max(nsumtj2)]
    wchr[i] <- names(nmark)[chind]
    chri <- ntj2[grep(wchr[i], bnams)]
    wint[i] <- (1:length(chri))[chri == max(chri)]
    nint <- names(chri)[wint[i]]
    qtl.fix[i] <- substr(nint, 6, nchar(nint))
    state[qtl.fix[i]] <- 0
    my.name <- paste("X",qtl.fix[i],sep="")
    message("Found QTL on chromosome ", wchr[i]," in interval ", wint[i])
    asdata[my.name] <- asdata[qtl.fix[i]]*100
    baseModel$call$data <- quote(asdata)
    add.qtl$call$data <- quote(asdata)
    baseModel$call$fixed <- as.formula(paste(deparse(baseModel$call$fixed, width.cutoff = 500), my.name, sep = " + "))
    baseModel$call$G.param <- baseModel$G.param
    baseModel$call$R.param <- baseModel$R.param    
    add.qtl$call$fixed <- baseModel$call$fixed
    add.qtl$call$G.param <- add.qtl$G.param
    add.qtl$call$R.param <- add.qtl$R.param
    i <- i + 1
    update <- TRUE
  }
  if(!is.null(wchr)) {
    baseModel$QTL$sig.chr <- substr(wchr, 5, nchar(wchr))
    baseModel$QTL$sig.int <- wint
    baseModel$QTL$diag$dmat <- dmat
    baseModel$QTL$diag$zrl <- zrl
    baseModel$QTL$diag$zcl <- zcl
    baseModel$QTL$diag$ochr <- ochr
    baseModel$QTL$diag$oint <- oint
}
  assign("asdata", asdata, envir = .GlobalEnv) 
  baseModel$QTL$qtlModel <- add.qtl
  class(baseModel) <- c("wgaim", "asreml")
  baseModel  
}

tr.wgaim <- function(object, iter = 1:length(object$QTL$sig.chr), diag.out = TRUE, ...){
     dots <- list(...)
     if(!is.na(pmatch("digits", names(dots))))
       dig <- dots$digits
     else
       dig <- options()$digits
     zrl <- object$QTL$diag$zrl
     if(any(ret <- is.na(pmatch(iter, 1:length(zrl))))){
      warning("\"iter\" values outside expected range .. using ones that are in iteration range")
      iter <- iter[!ret]
    }
     pvals <- unlist(lapply(zrl, function(el)
                           round(2*(1 - pnorm(abs(el))),dig)
                           ))       
     qtlmat <- matrix(NA, nrow = length(zrl), ncol = length(zrl))
     qtlmat[lower.tri(qtlmat, diag = TRUE)] <- pvals
     qnams <- paste(object$QTL$sig.chr, object$QTL$sig.int, sep = ".")                       
     dimnames(qtlmat) <- list(paste("Iter", 1:length(zrl), sep = "."), qnams)
     cat("\nIncremental QTL P-value Matrix.\n")
     cat("===============================\n")
     qtlmat[qtlmat < 0.001] <- "<0.001"
     qtlmat[is.na(qtlmat)] <- ""
     print.default(qtlmat[iter,1:iter[length(iter)]], quote = FALSE, right = TRUE, ...)
     if(diag.out) {
        cat("\nOutlier Detection Diagnostic.\n")
        cat("=============================\n")
        dmat <- round(as.matrix(object$QTL$diag$dmat), dig)
        dimnames(dmat)[[1]] <-  paste("Iter", 1:(length(zrl) + 1), sep = ".")
        dmat[,4][dmat[,4] < 0.001] <- "<0.001"
        print.default(dmat, quote = FALSE, right = TRUE, ...)
      }
  }

out.stat <- function (object, parentData, int = TRUE, iter = NULL, chr = NULL, 
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
  sc <- object$QTL$sig.chr
  si <- object$QTL$sig.int
  qtl <- paste(sc, si, sep = ".")
  dist <- lapply(parentData$geno, function(el) el$dist)
  ln <- length(sc)
  if (is.null(sc)) 
    stop("There are no outlier statistics to plot.")
  ochr <- object$QTL$diag$ochr
  oint <- object$QTL$diag$oint
  cnam <- names(ochr[[1]])
  cnam <- gsub("Chr.", "", cnam)
  inam <- names(oint[[1]])
  inam <- strsplit(gsub("qtls_Chr.", "", inam), "\\.")
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
  if (!int) 
    chr <- NULL
  if (!is.null(chr)) {
    dint <- dint[as.character(dint$chn) %in% chr, ]
    dist <- unlist(dist[chr])
  }
  else dist <- unlist(dist)
  labs <- unique(paste(dint$chn, ".", dint$inn, sep = ""))
  dist <- cumsum(dist) - dist/2
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
  if (int) {
    slabs <- strsplit(labs, "\\.")
    clabs <- unlist(lapply(slabs, function(el) el[1]))
    tabc <- table(clabs)
    cumc <- c(1, cumsum(tabc) + 1)
    labs[-cumc] <- " "
    labs <- as.character(labs)
    print(xyplot(int ~ dist | its, type = "l", data = dint, 
                 groups = chn, panel = panel.superpose,
                 panel.groups = function(x, y, ...) {
                   panel.xyplot(x, y, ...)
                 }, ylab = "outlier statistic", scales = list(x = list(labels = labs, 
                 rot = 45, at = unique(dist), cex = 0.6)), xlab = "Interval", 
                 layout = c(1, prow), ...))
  }
  else {
    dups <- !duplicated(paste(dint[, "chn"], dint[, c("its")], 
                              sep = "."))
    dchr <- dint[dups, ] 
    dchr$int <- unlist(ochr)
    names(dchr)[1] <- "chr"
    bc <- print(barchart(chr ~ chn | its, data = dchr, ylab = "outlier statistic", 
                   xlab = "Chromosome", layout = c(1, prow), ...))
    wchr <- dchr$chr[paste(dchr$chn, dchr$itn, sep = ":") %in% paste(td$chn, td$itn, sep = ":")]
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

print.wgaim <- function(x, parentData, ...){
  if(missing(parentData))
    stop("parentData is a required argument")
  if(!inherits(parentData, "cross"))
    stop("parentData is not of class \"cross\"")
  wchr <- x$QTL$sig.chr
  wint <- x$QTL$sig.int
  if(!length(x$QTL$sig.chr))
      cat("There are no significant putative QTL's\n")
  else {
    qtlm <- getQTL(x, parentData)
    for(z in 1:length(wchr)) {
      int <- paste(wchr[z], wint[z], sep = ".") 
      cat("\nPutative QTL found on the interval", int, "\nLeft-hand marker is", qtlm[z,1], "\nRight-hand marker is", qtlm[z,3], "\n")
    }
  }
}

## summary.wgaim <- function(object, parentData, ...){
##   if(missing(parentData))
##     stop("parentData is a required argument")
##   if(!inherits(parentData, "cross"))
##     stop("parentData is not of class \"cross\"")
##   if(is.null(object$QTL$sig.chr))
##     cat("There are no significant putative QTL's\n")
##   else {
##     wchr <- object$QTL$sig.chr
##     zrat <- round(object$coefficients$fixed/sqrt(object$vcoeff$fixed*object$sigma2), 2)
##     znam <- names(zrat)
##     zind <- asreml.grep("XChr\\.", znam)
##     zr <- rev(zrat[zind])
##     zc <- rev(object$coefficients$fixed[zind])
##     zn <- znam[zind]
##     zn <- substr(zn, 2, nchar(zn))
##     collab <- c("Chromosome", "Left Marker", "dist(cM)", "Right Marker", "dist(cM)", "Size", "z.ratio", "Pr(z)", "LOD")   
##     qtlmat <- matrix(ncol = 9, nrow = length(wchr))
##     qtlmat[,2:5] <- getQTL(object, parentData)
##     qtlmat[,1] <- wchr
##     qtlmat[,6:9] <- c(round(zc, 3), zr, round(2*(1 - pnorm(abs(zr))), 4),  round(0.5*log(exp(zr^2), base = 10), 4))
##     qtlmat <- matrix(qtlmat[order(qtlmat[,1]),], nrow = length(wchr))
##     dimnames(qtlmat) <- list(rowlab = as.character(1:length(wchr)), collab)
##     prmatrix(qtlmat, quote = FALSE, right = TRUE)
##     invisible(as.data.frame(qtlmat))
##   }
## }

summary.wgaim <- function(object, parentData, LOD = TRUE, ...){
  if(missing(parentData))
    stop("parentData is a required argument")
  if(!inherits(parentData, "cross"))
    stop("parentData is not of class \"cross\"")
  if(is.null(object$QTL$sig.chr))
    cat("There are no significant putative QTL's\n")
  else {
    wchr <- object$QTL$sig.chr
    zrat <- round(object$coefficients$fixed/sqrt(object$vcoeff$fixed*object$sigma2), 2)
    znam <- names(zrat)
    zind <- grep("X.", znam)
    zr <- rev(zrat[zind])
    zc <- rev(object$coefficients$fixed[zind])
    zn <- znam[zind]
    zn <- substr(zn, 2, nchar(zn))
    collab <- c("Chromosome", "Left Marker", "dist(cM)", "Right Marker", "dist(cM)", "Size", "z.ratio", "Pr(z)")   
    qtlmat <- matrix(ncol = 8, nrow = length(wchr))
    qtlmat[,2:5] <- getQTL(object, parentData)
    qtlmat[,1] <- wchr
    qtlmat[,6:8] <- c(round(zc, 3), zr, round(2*(1 - pnorm(abs(zr))), 4))
    if(LOD) {
      qtlmat <- cbind(qtlmat, round(0.5*log(exp(zr^2), base = 10), 4))
      collab <- c(collab, "LOD")
    }
    qtlmat <- matrix(qtlmat[order(qtlmat[,1]),], nrow = length(wchr))
    dimnames(qtlmat) <- list(as.character(1:length(wchr)), collab)
    class(qtlmat) <- "summary.wgaim"
    qtlmat
  }
}  

print.summary.wgaim <- function(x, ...){
  x <- unclass(x)
  print.default(x, quote = FALSE, right = TRUE, ...)
}
  
getQTL <- function(object, parentData){
  wchr <- object$QTL$sig.chr
  wint <- object$QTL$sig.int
  qtlm <- matrix(ncol = 4, nrow = length(wchr))
  for(i in 1:length(wchr)){
    rhmark <- parentData$geno[[wchr[i]]]$map[wint[i] + 1]
    lhmark <- parentData$geno[[wchr[i]]]$map[wint[i]]
    qtlm[i,1:4] <- c(names(lhmark),round(lhmark, 2),names(rhmark), round(rhmark, 2))
  }
 qtlm
}

## read.interval <- function(format = c("csv", "csvr", "csvs", "csvsr", "mm",
##    "qtx", "tlcart", "gary", "karl"), dir = getwd(), file, genfile, mapfile, 
##     phefile, chridfile, mnamesfile, pnamesfile, na.strings = c("-", 
##         "NA"), genotypes = c("A", "B"), estimate.map = TRUE, 
##     convertXdata = TRUE, missgeno = "MartinezCurnow",
##     rem.mark = TRUE, id = "id", subset = NULL, ...) {

##   oldops <- options()
##   options(warn = 1)
##   on.exit(options(oldops))
##   if(length(genotypes) > 2)
##     stop("This uses a modified version of the read.cross() function from
##   the qtl package. This strict modification only allows two genotypes \"A\" and \"B\"")
##   fullgeno <- read.crossQ(format = format, dir, file = file, genfile, mapfile, phefile, chridfile,
##    mnamesfile, pnamesfile, na.strings, genotypes = genotypes, convertXdata = convertXdata, estimate.map = estimate.map, ...)
##   fullgeno <- drop.nullmarkers(fullgeno)
##   if(!(id %in% names(fullgeno$pheno)))
##     stop("The unique identifier for the genotypic rows, ", deparse(substitute(id)), ",cannot be found in genotypic data") 
##   if(!is.null(subset))
##     fullgeno <- subset(fullgeno, chr = subset)
##   if(rem.mark) {
##     tpheno <- fullgeno$pheno
##     fullgeno <- fix.map(fullgeno, dir, file, id)
##     fullgeno$pheno <- tpheno
##   }
##   lid <- as.character(fullgeno$pheno[[id]])
##   fullgeno$geno <- lapply(fullgeno$geno, function(el, lid) {
##     row.names(el$data) <- as.character(lid)
##     el
##   }, lid)
##   chnam <- names(fullgeno$geno)
##   n.mark <- lapply(fullgeno$geno, function(x) length(x$map))
##   chnam.1 <- chnam[n.mark == 1]
##   chnam.2 <- chnam[n.mark != 1]
##   for (i in chnam.1) {
##         fullgeno$geno[[i]]$dist <- 0
##         fullgeno$geno[[i]]$theta <- 0
##         fullgeno$geno[[i]]$E.lambda <- 1/2
##         fullgeno$geno[[i]]$data[is.na(fullgeno$geno[[i]]$data)] <- 0
##         fullgeno$geno[[i]]$argmax <- fullgeno$geno[[i]]$data
##         for (j in 1:ncol(fullgeno$geno[[i]]$argmax)) {
##             fullgeno$geno[[i]]$argmax[, j] <- as.numeric(sub(2,
##                                                              -1, fullgeno$geno[[i]]$argmax[, j]))
##         }
##         fullgeno$geno[[i]]$intval <- as.matrix(fullgeno$geno[[i]]$argmax/2, ncol=1)
##         dimnames(fullgeno$geno[[i]]$intval)[[2]] <- names(fullgeno$geno[[i]]$map)
##       }
##   for(i in chnam.2){
##     fullgeno$geno[[i]]$dist <- diff(fullgeno$geno[[i]]$map)/100                                           
##     fullgeno$geno[[i]]$theta <- 0.5*(1-exp(-2*fullgeno$geno[[i]]$dist))                                        
##     fullgeno$geno[[i]]$E.lambda <- fullgeno$geno[[i]]$theta/(2*fullgeno$geno[[i]]$dist*(1-fullgeno$geno[[i]]$theta))                          
##   }
##   mtype <- c("Broman", "MartinezCurnow")
##   if(is.na(type <- pmatch(missgeno, mtype)))
##     stop("Missing marker type must be one of \"Broman\" or \"MartinezCurnow\". Partial matching is allowed.")
##   missgeno <- mtype[type]  
##   if(missgeno=="MartinezCurnow"){
##     for(i in chnam.2){
##       dtemp <- fullgeno$geno[[i]]$data
##       for(j in 1:ncol(fullgeno$geno[[i]]$data))
##         dtemp[,j] <-  as.numeric(sub(2,-1,fullgeno$geno[[i]]$data[,j]))
##       fullgeno$geno[[i]]$argmax <- dtemp
##     }
##     for(p in 1:length(chnam.2))
##       fullgeno$geno[[p]]$argmax <- miss.q(fullgeno$geno[[p]]$theta, fullgeno$geno[[p]]$argmax)
##   }
##   else {
##     fullgeno <- argmax.geno(fullgeno)
##     for(i in chnam.2){
##       for(j in 1:ncol(fullgeno$geno[[i]]$argmax))
##         fullgeno$geno[[i]]$argmax[,j] <- as.numeric(sub(2,-1,fullgeno$geno[[i]]$argmax[,j]))
##       dimnames(fullgeno$geno[[i]]$argmax)[[1]] <- dimnames(fullgeno$geno[[1]]$data)[[1]]
##     }    
##   }
##   for(i in chnam.2){
##     lambda <- addiag(fullgeno$geno[[i]]$E.lambda,-1) + addiag(c(fullgeno$geno[[i]]$E.lambda,0),0) 
##     lambda <- lambda[,-dim(lambda)[2]]
##     fullgeno$geno[[i]]$intval <- fullgeno$geno[[i]]$argmax %*% lambda
##     dimnames(fullgeno$geno[[i]]$intval)[[2]] <- names(fullgeno$geno[[i]]$theta)
##   }         
##   class(fullgeno) <- c(class(fullgeno), "interval")
##   fullgeno
## }

read.interval <- function(format = c("csv", "csvr", "csvs", "csvsr", "mm",
   "qtx", "tlcart", "gary", "karl"), dir = getwd(), file, genfile, mapfile, 
    phefile, chridfile, mnamesfile, pnamesfile, na.strings = c("-", 
        "NA"), genotypes = c("A", "B"), estimate.map = TRUE, 
    convertXdata = TRUE, missgeno = "MartinezCurnow",
    rem.mark = TRUE, id = "id", subset = NULL, ...) {
  cat("Users should be aware that this function will be deprecated in future releases of \"wgaim\".\n") 
  oldops <- options()
  options(warn = 1)
  on.exit(options(oldops))
  if(length(genotypes) > 2)
    stop("This uses a modified version of the read.cross() function from
  the qtl package. This strict modification only allows two genotypes, for example, \"A\" and \"B\"")
  fullgeno <- read.crossQ(format = format, dir, file = file, genfile, mapfile, phefile, chridfile,
   mnamesfile, pnamesfile, na.strings, genotypes = genotypes, convertXdata = convertXdata, estimate.map = estimate.map, ...)
  cross2int(fullgeno = fullgeno, missgeno = missgeno, rem.mark = rem.mark, id = id, subset = subset)
 } 

cross2int <- function(fullgeno, missgeno = "MartinezCurnow", rem.mark = TRUE, id = "id",
                      subset = NULL){
  if(!inherits(fullgeno, "bc"))
    stop("This function is restricted to population containing only two genotypes.")
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
  fullgeno$geno <- lapply(fullgeno$geno, function(el, lid) {
    row.names(el$data) <- as.character(lid)
    el
  }, lid)
  chnam <- names(fullgeno$geno)
  n.mark <- lapply(fullgeno$geno, function(x) length(x$map))
  chnam.1 <- chnam[n.mark == 1]
  chnam.2 <- chnam[n.mark != 1]
  for (i in chnam.1) {
        fullgeno$geno[[i]]$dist <- 0
        fullgeno$geno[[i]]$theta <- 0
        fullgeno$geno[[i]]$E.lambda <- 1/2
        fullgeno$geno[[i]]$data[is.na(fullgeno$geno[[i]]$data)] <- 0
        fullgeno$geno[[i]]$argmax <- fullgeno$geno[[i]]$data
        for (j in 1:ncol(fullgeno$geno[[i]]$argmax)) {
            fullgeno$geno[[i]]$argmax[, j] <- as.numeric(sub(2,
                                                             -1, fullgeno$geno[[i]]$argmax[, j]))
        }
        fullgeno$geno[[i]]$intval <- as.matrix(fullgeno$geno[[i]]$argmax/2, ncol=1)
        dimnames(fullgeno$geno[[i]]$intval)[[2]] <- names(fullgeno$geno[[i]]$map)
      }
  for(i in chnam.2){
    fullgeno$geno[[i]]$dist <- diff(fullgeno$geno[[i]]$map)/100                                           
    fullgeno$geno[[i]]$theta <- 0.5*(1-exp(-2*fullgeno$geno[[i]]$dist))                                        
    fullgeno$geno[[i]]$E.lambda <- fullgeno$geno[[i]]$theta/(2*fullgeno$geno[[i]]$dist*(1-fullgeno$geno[[i]]$theta))                          
  }
  mtype <- c("Broman", "MartinezCurnow")
  if(is.na(type <- pmatch(missgeno, mtype)))
    stop("Missing marker type must be one of \"Broman\" or \"MartinezCurnow\". Partial matching is allowed.")
  missgeno <- mtype[type]  
  if(missgeno=="MartinezCurnow"){
    for(i in chnam.2){
      dtemp <- fullgeno$geno[[i]]$data
      for(j in 1:ncol(fullgeno$geno[[i]]$data))
        dtemp[,j] <-  as.numeric(sub(2,-1,fullgeno$geno[[i]]$data[,j]))
      fullgeno$geno[[i]]$argmax <- dtemp
    }
    for(p in 1:length(chnam.2))
      fullgeno$geno[[p]]$argmax <- miss.q(fullgeno$geno[[p]]$theta, fullgeno$geno[[p]]$argmax)
  }
  else {
    fullgeno <- argmax.geno(fullgeno)
    for(i in chnam.2){
      for(j in 1:ncol(fullgeno$geno[[i]]$argmax))
        fullgeno$geno[[i]]$argmax[,j] <- as.numeric(sub(2,-1,fullgeno$geno[[i]]$argmax[,j]))
      dimnames(fullgeno$geno[[i]]$argmax)[[1]] <- dimnames(fullgeno$geno[[1]]$data)[[1]]
    }    
  }
  for(i in chnam.2){
    lambda <- addiag(fullgeno$geno[[i]]$E.lambda,-1) + addiag(c(fullgeno$geno[[i]]$E.lambda,0),0) 
    lambda <- lambda[,-dim(lambda)[2]]
    fullgeno$geno[[i]]$intval <- fullgeno$geno[[i]]$argmax %*% lambda
    dimnames(fullgeno$geno[[i]]$intval)[[2]] <- names(fullgeno$geno[[i]]$theta)
  }         
  class(fullgeno) <- c(class(fullgeno), "interval")
  fullgeno
}

wmerge <- function (geno, pheno, by = NULL, ...)
{
    if (missing(pheno) | missing(geno))
        stop("geno and pheno are required arguments")
    if (!inherits(geno, "interval"))
        stop("Genotypic data needs to be of class \"interval\".")
    if (is.null(by))
        stop("Need name of matching column to merge datasets.")
    if (is.character(pheno))
        pheno <- asreml.read.table(pheno, ...)
    if (is.null(other <- geno$pheno[, by]))
        stop("Genotypic data does not contain column \"", by, "\".")
    if (is.null(pheno[, by]))
        stop("Phenotypic data does not contain column \"", by, "\".")
    mby <- pmatch(as.character(other), as.character(pheno[, by]))
    if (all(is.na(mby)))
        stop("Names in Genotypic \"", by, "\" column do not match any names in Phenotypic \"",
            by, "\" column.")
    gdat <- lapply(geno$geno, function(el) el$intval)
    n.int <- lapply(geno$geno, function(el) dim(el$intval)[2])
    temp <- cbind.data.frame(other, do.call("cbind", gdat))
    nch <- rep(names(nmar(geno)), unlist(n.int))
    nint <- list()
    for (i in 1:nchr(geno)) nint[[i]] <- seq(1, n.int[[i]])
    nint <- unlist(nint)
    names(temp) <- c(by, paste("Chr", nch, nint, sep = "."))
    geno$pheno.dat <- pheno
    pheno <- cbind(ord = 1:nrow(pheno), pheno)
    geno$full.data <- merge(pheno, temp, by.x = by, by.y = by,
        all.x = TRUE, all.y = FALSE)
    geno$full.data <- geno$full.data[order(geno$full.data$ord),]
    geno$full.data <- geno$full.data[, -2]
    geno
}


read.crossQ <- function (format = c("csv", "csvr", "csvs", "csvsr", "mm",
   "qtx", "tlcart", "gary", "karl"), dir, file, genfile, mapfile, 
    phefile, chridfile, mnamesfile, pnamesfile, na.strings, genotypes,
    estimate.map, convertXdata, ...) 
{
    if (format == "csvrs") {
        format <- "csvsr"
        warning("Assuming you mean 'csvsr' rather than 'csvrs'.\n")
    }
    format <- match.arg(format)
    if (format == "csv" || format == "csvr") {
        cross <- read.cross.csv(dir, file, na.strings, genotypes, 
            estimate.map, rotate = (format == "csvr"), ...)
    }
    else if (format == "csvs" || format == "csvsr") {
        if (missing(phefile) && !missing(file) && !missing(genfile)) {
            phefile <- genfile
            genfile <- file
        }
        else if (missing(genfile) && !missing(file) && !missing(phefile)) {
            genfile <- file
        }
        cross <- read.cross.csvs(dir, genfile, phefile, na.strings, 
            genotypes, estimate.map, rotate = (format == "csvsr"), 
            ...)
    }
    else if (format == "qtx") {
        cross <- read.cross.qtx(dir, file, estimate.map)
    }
    else if (format == "qtlcart") {
        if (missing(mapfile) && !missing(genfile)) 
            mapfile <- genfile
        cross <- read.cross.qtlcart(dir, file, mapfile)
    }
    else if (format == "karl") {
        if (missing(genfile)) 
            genfile <- "gen.txt"
        if (missing(mapfile)) 
            mapfile <- "map.txt"
        if (missing(phefile)) 
            phefile <- "phe.txt"
        cross <- read.cross.karl(dir, genfile, mapfile, phefile)
    }
    else if (format == "mm") {
        if (missing(mapfile) && !missing(genfile)) 
            mapfile <- genfile
        cross <- read.cross.mm(dir, file, mapfile, estimate.map)
    }
    else if (format == "gary") {
        if (missing(genfile)) 
            genfile <- "geno.dat"
        if (missing(mnamesfile)) 
            mnamesfile <- "mnames.txt"
        if (missing(chridfile)) 
            chridfile <- "chrid.dat"
        if (missing(phefile)) 
            phefile <- "pheno.dat"
        if (missing(pnamesfile)) 
            pnamesfile <- "pnames.txt"
        if (missing(mapfile)) 
            mapfile <- "markerpos.txt"
        cross <- read.cross.gary(dir, genfile, mnamesfile, chridfile, 
            phefile, pnamesfile, mapfile, estimate.map, na.strings)
    }
    estimate.map <- cross[[2]]
    cross <- cross[[1]]
    chrnam <- names(cross$geno)
    if (all(regexpr("^[Cc][Hh][Rr]", chrnam) > 0)) {
        chrnam <- substr(chrnam, 4, nchar(chrnam))
        if (all(regexpr("^[Oo][Mm][Oo][Ss][Oo][Mm][Ee]", chrnam) > 
            0)) 
            chrnam <- substr(chrnam, 8, nchar(chrnam))
    }
    if (sum(chrnam == "x") > 0) 
        chrnam[chrnam == "x"] <- "X"
    names(cross$geno) <- chrnam
    for (i in 1:length(cross$geno)) if (names(cross$geno)[i] == 
        "X") 
        class(cross$geno[[i]]) <- "X"
    chrtype <- sapply(cross$geno, class)
    if (any(chrtype == "X") && convertXdata) {
        if (class(cross)[1] == "bc") 
            cross <- fixXgeno.bc(cross)
        if (class(cross)[1] == "f2") 
            cross <- fixXgeno.f2(cross)
    }
    if (estimate.map) {
        cat(" --Estimating genetic map\n")
        newmap <- est.map(cross)
        cross <- replace.map(cross, newmap)
    }
    for (i in 1:nchr(cross)) storage.mode(cross$geno[[i]]$data) <- "integer"
    summary(cross)
    cross
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
      if((n %in% whk)) {
        if(length(whk) >= 2)
          whk <- c(whk[length(whk)], whk[-length(whk)])
      }
#        { 
#        if(length(whk) > 2)
#          whk <- c(whk[1], whk[length(whk)], whk[2:(length(whk) - 1)])
#        if(length(whk) == 2)
#          whk <- c(whk[length(whk)], whk[1])
#      }
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
    chr
  } 

link.map <- function(object, ...)
  UseMethod("link.map")

link.map.cross <- function(object, chr, max.dist, marker.names = "markers", tick = FALSE, squash = TRUE, ...){

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
      cex <- dots$cex
    else
      cex <- par("cex")
    if(!squash)
      chrpos <- seq(1, n.chr * 3, by = 3)
    else
      chrpos <- seq(1, n.chr * 2, by = 2)
    thelim <- range(chrpos) + c(-1.6, 1.35)
    for(i in 1:n.chr) {
      mt[[i]] <- map[[i]]
      conv <- par("pin")[2]/maxlen    
      for(j in 1:(length(mt[[i]]) - 1)){
        ch <- mt[[i]][j + 1]*conv - (mt[[i]][j]*conv + 10*par("csi")*cex/9)
        if(ch < 0){
          temp <- mt[[i]][j + 1]*conv + abs(ch)
          mt[[i]][j + 1] <- temp/conv       
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
         text(chrpos[i] + 0.50, mt[[i]], as.character(round(map[[i]], 2)), 
              adj = c(0, 0.5), ...)
       else 
         text(chrpos[i] + 0.50, mt[[i]], names(map[[i]]), adj = c(0, 0.5), ...)
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

circ <- function(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
    ((x - shiftx)^2)/elx + ((y - shifty)^2)/ely

link.map.default <- function(object, parentData, chr, max.dist, marker.names = "markers", clist = list(qcol = rainbow(length(object)), mcol = "red", tcol = rainbow(length(object))), trait.labels = NULL, tick = FALSE, ...){

  old.par <- par(no.readonly = TRUE)
  par(mar = c(5, 4, 5, 2) + 0.1)
  par(xpd = TRUE)
  on.exit(par(old.par))
  dlist <- list()
  lclass <- lapply(object, function(el){
    if(!inherits(el, "wgaim"))
      stop("list objects need to inherit from class \"wgaim\"")
    class(el)[1]
  })
  lclass <- unique(unlist(lclass))
  if(any(is.na(pmatch(lclass, "wgaim"))))
      stop("link.map method is only for list objects of class \"wgaim\"")
  sig.chr <- lapply(object, function(el) el$QTL$sig.chr)
  len <- lapply(sig.chr, length)
  if(!is.null(trait.labels)) {
    if(length(trait.labels) != length(object))
      stop("Length of trait labels does not equal number of models specified.")
    trait <- rep(trait.labels, times = len)
  }
  else { 
    trait <- unlist(lapply(object, function(el)
        rep(as.character(el$call$fixed[[2]]), length(el$QTL$sig.chr))))
  } 
  dlist$QTL$sig.chr <- unlist(sig.chr) 
  dlist$QTL$sig.int <- unlist(lapply(object, function(el) el$QTL$sig.int))
#  attr(dlist, "trait") <- trait                       
  class(dlist) <- "wgaim"
  if(length(clist$qcol) != length(object)) {
    warning("QTL colours not of the same length as the number of traits,\n Choosing \"qcol = rainbow(",length(object),")\".")
    clist$qcol <- rainbow(length(object))
  }
  if(length(clist$tcol) != length(object)){
    if(length(clist$tcol) != 1) {
      warning("Inappropriate length for trait name colours, using QTL colours")
      clist$tcol <- clist$qcol
    }
  }
  link.map(dlist, parentData, chr, max.dist, marker.names = marker.names, clist = clist, trait.labels = trait, tick = tick, ...)
#  legend(x = par("usr")[2]/2, y = par("usr")[4]*2.2, unique(trait), fill = qcol, xjust = 0.5, horiz = TRUE, bty = "n")
 }

link.map.wgaim <- function(object, parentData, chr, max.dist, marker.names = "markers", clist = list(qcol = "light blue", mcol = "red", tcol = "light blue"), trait.labels = NULL, tick = FALSE, ...){   

  dots <- list(...)
  if(missing(parentData))
    stop("parentData is a required argument")
  if(!inherits(parentData, "cross"))
    stop("parentData is not of class \"cross\"")
  if(!length(wchr <- object$QTL$sig.chr)){
    warning("There are no significant QTL's. Plotting map only...")
    link.map(parentData, chr, max.dist, marker.names = marker.names, tick = tick, squash = FALSE, ...)
    return(invisible())
  }
  if(is.null(clist$qcol))
    clist$qcol <- "light blue"
  if(is.null(clist$tcol))
    clist$tcol <- clist$qcol
  if(missing(chr))
    chr <- unique(wchr)[order(unique(wchr))]
  lmap <- link.map(parentData, chr, max.dist, marker.names = marker.names, tick = tick, squash = TRUE, ...)
  map <- lmap$map
  qtlm <- getQTL(object, parentData)
  if(is.null(trait <- trait.labels))
    trait <- rep(as.character(object$call$fixed[[2]]), length(wchr))
  if(length(trait.labels) == 1)
    trait <- rep(trait.labels, length(wchr))
  qtrait <- unique(trait) 
  qtlm <- cbind.data.frame(qtlm, trait = factor(trait, levels = unique(trait)))
  qtlm[,2] <- as.numeric(as.character(qtlm[,2]))
  qtlm[,4] <- as.numeric(as.character(qtlm[,4])) 
  qtlList <- lapply(split(qtlm, wchr), function(el) el[order(el[,2]),])
  qtlm <- do.call("rbind", qtlList)
  wchr <- wchr[order(wchr)]
  if(!missing(chr)) {
    if(any(is.na(wh <- pmatch(wchr, chr, dup = TRUE)))){
      warning("Some QTL's exist outside chromosome(s) subset, Omitting QTL's....")
      qtlm <- qtlm[!is.na(wh),]
      wchr <- wchr[!is.na(wh)]
#      trait <- trait[!is.na(wh)]
    }  
  }
  if(!missing(max.dist)) {
    rml  <- qtlm[,4] > max.dist
    if(any(rml)){
      warning("Some QTL regions outside maximum distance specified. Omitting QTL's....") 
#      qtlm <- matrix(qtlm[!rml,], nrow = length(rml[!rml]), byrow = FALSE)
      qtlm <- qtlm[!rml,]
      wchr <- wchr[!rml]
#      trait <- trait[!rml]
    }
  }
  n.chr <- length(map)
  maxlen <- max(unlist(lapply(map, max)))
  chrpos <- lmap$chrpos
  mt <- lmap$mt
  if(!is.na(cind <- pmatch("col", names(dots))))
    dots <- dots[-cind]
  if(is.null(dim(qtlm)))
    qtlm <- matrix(qtlm, nrow = 1, byrow = FALSE)

  ## move trait names on map for multiple qtls

  tlis <- list()
  if(!is.na(pmatch("cex", names(dots))))
    cex <- dots$cex
  else
    cex <- par("cex")
  for(i in 1:n.chr){
    conv <- par("pin")[2]/maxlen        
    if(as.logical(length(ind <- grep(names(map)[i], wchr)))){
      tlis[[i]] <- as.vector(qtlm[ind, 2] + qtlm[ind,4])/2
      names(tlis[[i]]) <- as.character(qtlm[ind,5])
      if(length(ind) > 1) {          
          for(j in 1:(length(ind) - 1)){
            ch <- tlis[[i]][j + 1]*conv - (tlis[[i]][j]*conv + 10*par("csi")*cex/9)
            if(ch < 0){
              temp <- tlis[[i]][j + 1]*conv + abs(ch)
              tlis[[i]][j + 1] <- temp/conv       
            }
          }
        }
    }
  }
  ## check duplication in marker names for multiple traits 
  
  qtld <- qtlm[,1:4]
  nodup <- !duplicated(do.call("paste",qtld))
  qtls <- qtld[nodup,]
  whd <- pmatch(do.call("paste", qtld), do.call("paste", qtls), dup = TRUE)
  dlis <- split(as.character(qtlm[,5]), whd)
  qtlm <- qtls
  wchr <- wchr[nodup]
  for(i in 1:n.chr){
    if(as.logical(length(ind <- grep(names(map)[i], wchr)))){
       for(j in ind) {
        if(!is.null(marker.names)){
          wh <- mt[[i]][pmatch(c(as.character(qtlm[j,1]), as.character(qtlm[j,3])), names(map[[i]]))]
          if(marker.names == "dist"){
            dist <- map[[i]][pmatch(c(as.character(qtlm[j,1]), as.character(qtlm[j,3])), names(map[[i]]))]
            alis <- list(x = chrpos[i] + 0.50, y = wh, labels = as.character(round(as.numeric(dist), 2)), adj = c(0, 0.5), col = clist$mcol)
          }
          else
            alis <- list(x = chrpos[i] + 0.50, y = wh, labels = names(wh), adj = c(0, 0.5), col = clist$mcol)
          do.call("text", c(alis, dots))
        }
        yv <- c(qtlm[j,2],qtlm[j,4])
        yv <- c(yv, rev(yv))
        dind <- dlis[[j]]
        qcols <- clist$qcol[pmatch(dind, qtrait)]
        qind <- 1:length(dind)
        if(length(dlis[[j]]) > 1) {
            int <- seq(chrpos[i] - 0.20, chrpos[i] + 0.20, length = length(dind) + 1)
            for(k in 1:length(dind)){ 
              xv <- c(rep(int[k], 2), rep(int[k + 1], 2)) 
              polygon(xv, y = yv, border = NA, col = qcols[k])
          }                   
        }
        else {
          xv <- c(rep(chrpos[i] - 0.20, 2), rep(chrpos[i] + 0.20, 2))       
          polygon(xv, y = yv, border = NA, col = qcols)
      }
        segments(chrpos[i] - 0.25, yv[1], chrpos[i] - 0.25, yv[2])
        segments(chrpos[i] - 0.25, sum(yv[1:2])/2, chrpos[i] - 0.3, sum(yv[1:2])/2)        
        segments(chrpos[i] - 0.3, sum(yv[1:2])/2, chrpos[i] - 0.4,  tlis[[i]][qind])
        segments(chrpos[i] - 0.4, tlis[[i]][qind], chrpos[i] - 0.45, tlis[[i]][qind])
        if(length(clist$tcol) > 1)
          tcols <- clist$tcol[pmatch(dind, qtrait)]
        else
          tcols <- clist$tcol
        text(chrpos[i] - 0.5, tlis[[i]][qind], names(tlis[[i]][qind]), adj = c(1, 0.3), col = tcols, cex = cex)
        tlis[[i]] <- tlis[[i]][-qind]
      }
    }
    segments(chrpos[i] - 0.20, map[[i]], chrpos[i] + 
             0.20, map[[i]])
  }
  if(is.na(pmatch("main", names(dots))))
     title("Genetic Map with QTLs")
}







