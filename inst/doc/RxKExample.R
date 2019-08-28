## ---- include = FALSE----------------------------------------------------
library(wgaim)
library(ggplot2)
knitr::opts_chunk$set(fig.path = "", highlight = FALSE)

## ------------------------------------------------------------------------
data(package = "wgaim")

## ------------------------------------------------------------------------
wgpath <- system.file("extdata", package = "wgaim")
list.files(wgpath)

## ------------------------------------------------------------------------
data(phenoRxK, package = "wgaim")
head(phenoRxK)

## ------------------------------------------------------------------------
rkyld.asi <- asreml::asreml(yld ~ Type, random = ~ Genotype + Rep, residual =
     ~ ar1(Range):ar1(Row), data = phenoRxK)

## ------------------------------------------------------------------------
summary(rkyld.asi)$varcomp

## ----variogram, fig.align = "center", out.width = "95%", dpi = 150-------
plot(asreml::varioGram.asreml(rkyld.asi))

## ----residuals, warning = FALSE, fig.align = "center", fig.width = 13, fig.height = 10, out.width = "95%", dpi = 150----
phenoRxKd <- cbind.data.frame(phenoRxK, Residuals = resid(rkyld.asi))
ggplot(phenoRxKd, aes(y = Residuals, x = as.numeric(Range))) + facet_wrap(~ Row) +
    geom_hline(yintercept = 0, linetype = 2) + geom_point(shape = 16, colour = "blue") +
    xlab("Range") + theme_bw()

## ------------------------------------------------------------------------
rkyld.asf <- asreml::asreml(yld ~ Type + lrow, random = ~ Genotype + Range, residual =
     ~ ar1(Range):ar1(Row), data = phenoRxK)

## ------------------------------------------------------------------------
data(genoRxK, package = "wgaim")

## ------------------------------------------------------------------------
genoRxK <- read.cross(format = "csvr", file="genoRxK.csv", genotypes=c("AA","BB"),
     dir = wgpath, na.strings = c("-", "NA"))

## ------------------------------------------------------------------------
summary(genoRxK)
names(genoRxK$pheno)

## ------------------------------------------------------------------------
genoRxKi <- cross2int(genoRxK, consensus.mark = TRUE, impute = "MartinezCurnow",
                      id = "Genotype")

## ------------------------------------------------------------------------
names(genoRxKi$geno[[1]])

## ---- warning = FALSE, eval = FALSE--------------------------------------
#  rkyld.qtl <- wgaim(rkyld.asf, genoRxKi, merge.by = "Genotype", fix.lines = TRUE,
#                      gen.type = "interval", method = "random", selection = "interval",
#                      trace = "rxk.txt", na.action = asreml::na.method(x = "include"))

## ---- warning = FALSE, eval = TRUE---------------------------------------
rkyld.asf <- asreml::update.asreml(rkyld.asf, fixed. = . ~ . - Type)
rkyld.qtl <- wgaim(rkyld.asf, genoRxKi, merge.by = "Genotype", fix.lines = TRUE,
                    gen.type = "interval", method = "random", selection = "interval",
                    trace = "rxk.txt", na.action = asreml::na.method(x = "include"))

## ------------------------------------------------------------------------
names(rkyld.qtl$QTL)
names(rkyld.qtl$QTL$diag)
class(rkyld.qtl)

## ------------------------------------------------------------------------
summary(rkyld.qtl, genoRxKi)

## ----out, warning = FALSE, fig.align = "center", fig.width = 13, fig.height = 6, out.width = "95%", dpi = 150----
outStat(rkyld.qtl, genoRxKi, iter = 1:2, statistic = "outlier")
outStat(rkyld.qtl, genoRxKi, iter = 1:2, statistic = "blups")

## ----link, fig.align = "center", fig.width = 13, fig.height = 6, out.width = "95%", dpi = 150----
linkMap(rkyld.qtl, genoRxKi)

