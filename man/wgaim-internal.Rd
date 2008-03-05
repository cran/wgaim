\name{wgaim-internal}
\alias{fix.map}
\alias{circ}
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
circ(x, y, shiftx = 0, shifty = 0, ely = 1, elx = 1)
getQTL(object, parentData)
addiag(x = 1, di = 0, nrow.arg, ncol.arg = n)
miss.q(theta, chr)
read.crossQ(format = c("csv", "csvr", "csvs", "csvsr", "mm",
   "qtx", "tlcart", "gary", "karl"), dir, file, genfile, mapfile,
phefile, chridfile, mnamesfile, pnamesfile, na.strings, genotypes,
estimate.map, convertXdata, \ldots)

\method{wgaim}{default}(baseModel, \ldots) 
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\examples{}
\keyword{internal}
