\name{wgaim-internal}
\alias{fix.map}
\alias{link.map.default}
\alias{getQTL}
\alias{addiag}
\alias{miss.q}
\alias{read.crossQ}
\alias{wgaim.default}
\title{Internal wgaim functions}
\description{
These are internal \code{wgaim} functions that are not intended to be directly called by the user  
}
\usage{
fix.map(full.data, dir, file, id)
getQTL(object, parentData)
addiag(x = 1, di = 0, nrow.arg, ncol.arg = n)
miss.q(theta, chr)
read.crossQ(format = c("csv", "csvr", "csvs", "csvsr", "mm",
   "qtx", "tlcart", "gary", "karl"), dir, file, genfile, mapfile,
phefile, chridfile, mnamesfile, pnamesfile, na.strings, genotypes,
estimate.map, convertXdata, \ldots)
\method{link.map}{default}(object, parentData, chr, max.dist, marker.names = TRUE,
qcol = rainbow(length(object)), mcol = "red", trait.labels = NULL, \ldots)
\method{wgaim}{default}(baseModel, \ldots) 
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\examples{}
\keyword{internal}
