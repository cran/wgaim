\name{wgaim-internal}
\alias{fix.map}
\alias{getQTL}
\alias{addiag}
\alias{miss.q}
\alias{wgaim.default}
\alias{qtl.pick}
\alias{mergeData}
\alias{updateWgaim}
\alias{error.code}
\alias{v.modify}
\alias{v.init}
\title{Internal wgaim functions}
\description{
These are internal \code{wgaim} functions that are not intended to be directly called by the user
}
\usage{
fix.map(full.data, dir, filename, id, \ldots)
getQTL(object, intervalObj)
addiag(x = 1, di = 0, nrow.arg, ncol.arg = n)
miss.q(theta, chr)
qtl.pick(asr, intervalObj, gen.type, selection, state, verboseLev)
mergeData(phenoData, geneticData, by)
error.code(ec = NULL)
v.modify(model)
v.init(term, G.param)
updateWgaim(object, attempts, \ldots)
\method{wgaim}{default}(baseModel, \ldots)
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\keyword{internal}
