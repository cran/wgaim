## 2012 wgaim 1.3+ script

########### Sunco x Tasman Milling Experiment
########### Example 2:

## import phenotypic data

wgpath <- system.file("extdata", package = "wgaim")
phenoSxT <- asreml.read.table(paste(wgpath,"\\phenoSxT.csv", sep =""), header = TRUE, sep = ",")

## initial model and diagnostic

st.fmI <- asreml(myield ~ Type, random = ~ id + Rep + Range:Row + Millday, rcov = ~ Millday:ar1(Millord), data = phenoSxT)
xyplot(resid(st.fmI) ~ as.numeric(Millord) | Millday, data = phenoSxT, type = "b", xlab = "Milling Order", ylab = "Residuals")
field.resid <- coef(st.fmI, pattern = "Range:Row")
rrd <- data.frame(field.resid = field.resid, Range = factor(rep(1:12, each = 31)), Row = rep(1:31, 12))
xyplot(field.resid ~ Row | Range, fun = "xyplot", data = rrd, type = "b")
summary(st.fmI)$varcomp

# full base model

st.fmF <- asreml(myield ~ Type + lord + lrow, random = ~ id + Rep + Range:Row + Millday, rcov = ~ Millday:ar1(Millord), data = phenoSxT, na.method.X = "include")
summary(st.fmF)$varcomp

# fit null model (for comparison)

st.fmN <- asreml(myield ~ 1, random = ~ id, data = phenoSxT, na.method.X = "include")

## read in genetic data and check

genoSxT <- read.cross("csvr", file="genoSxT.csv", genotypes=c("AA","BB"),
                    dir = wgpath, na.strings = c("-", "NA"))
nmar(genoSxT)
nt <- ntyped(genoSxT, "ind")
nt[nt < 120]
genoSxT <- subset(genoSxT, ind = 1:180)
genoSxT <- cross2int(genoSxT, missgeno="Mart", id = "id", rem.mark = TRUE, map.function = "kosambi")
nmar(genoSxT)

# plot the linkage map in various ways

link.map(genoSxT, chr = names(nmar(genoSxT)), m.cex = 0.5)
link.map(genoSxT, names(nmar(genoSxT)), m.cex = 0.5,
         chr.dist = list(start = 25, end = 180), marker.names = "dist")

## wgaim QTL analyses for Full and Null model

st.qtlN <- wgaim(st.fmN, phenoData = phenoSxT, intervalObj = genoSxT,
    merge.by = "id", gen.type = "interval", method = "fixed",
    selection = "interval", trace = "nullmodel.txt")

st.qtlF <- wgaim(st.fmF, phenoData = phenoSxT, intervalObj = genoSxT,
    merge.by = "id", gen.type = "interval", method = "fixed",
    selection = "interval", trace = "fullmodel.txt")

## diagnostics

# plot the outlier statistics

out.stat(st.qtlF, genoSxT, int = TRUE, iter = 1:2, cex = 0.6, stat = "os")
out.stat(st.qtlF, genoSxT, int = TRUE, iter = 1:2, cex = 0.6, stat = "blups")
out.stat(st.qtlF, genoSxT, int = TRUE, iter = 1:5, cex = 0.6, chr = c("2B","4B","6B","7D"))

# trace of forward selection proces

tr(st.qtlF, iter = 1:10, digits = 3)

## visualisation and summary

# summary and table

summary(st.qtlF, genoSxT, LOD = TRUE)
qtlTable(st.qtlF, st.qtlN, intervalObj = genoSxT, labels = c("Full", "Null"), columns = 1:8)

# plot qtl on linkage map in various ways

link.map(st.qtlF, genoSxT, marker.names = "dist", trait.labels = "Full")
link.map(st.qtlF, genoSxT, chr = c("1B1", "2B", "5D"), marker.names = "dist", trait.labels = "Full")
link.map(st.qtlF, genoSxT, chr.dist = list(start = 25), marker.names = "dist", trait.labels = "Full")

# with both models

link.map.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
                 trait.labels = c("Full", "Null"), list.cex = list(m.cex = 0.7))
link.map.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
                 trait.labels = c("Full", "Null"), chr.dist = list(start = 25))

# customize your qtl plots

link.map.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = c("skyblue3",
    "salmon2"), m.col = "red", t.col = c("skyblue3", "salmon2")), col = "gray")

link.map.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = rep(gray(0.8), 2),
    m.col = "black", t.col = "black"), list.cex = list(t.cex = 0.8), col = "gray")

## whole genome marker analysis and summary

st.qtlFM <- wgaim(st.fmF, phenoData = phenoSxT, intervalObj = genoSxT,
                 merge.by = "id", gen.type = "marker", method = "fixed",
                 selection = "interval", trace = "fullmodel.txt")

summary(st.qtlFM, genoSxT, LOD = TRUE)
out.stat(st.qtlFM, genoSxT, int = TRUE, iter = 1:5, cex = 0.6)
out.stat(st.qtlFM, genoSxT, int = TRUE, iter = 1:5, cex = 0.6, chr = c("2B","4B","5A", "6B","7D"))
link.map(st.qtlFM, genoSxT, marker.names = "dist", trait.labels = "Full", cex = 2.0, pch = 16)

# Null model

st.qtlNM <- wgaim(st.fmN, phenoData = phenoSxT, intervalObj = genoSxT,
                 merge.by = "id", gen.type = "marker", method = "fixed",
                 selection = "interval", trace = "nullmodel.txt")

# linkage map with QTL

link.map.default(list(st.qtlFM, st.qtlNM), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = c("red",
    "light blue"), m.col = "red", t.col = c("red", "light blue")),
    list.cex = list(t.cex = 0.9, m.cex = 0.7), col = "black", cex = 2, pch = 16)

link.map.default(list(st.qtlFM, st.qtlNM), genoSxT, marker.names = "dist",
    chr.dist = list(start = 50), trait.labels = c("Full", "Null"),
    list.col = list(q.col = c("red", "light blue"), m.col = "red", t.col
    = c("red", "light blue")), list.cex = list(t.cex = 0.9, m.cex = 0.7),
    col = "black", cex = 2, pch = 16)


########### end script
