\name{wgaim-package}
\alias{wgaim-package}
\docType{package}
\title{Whole Genome Average Interval Mapping (wgaim) for QTL detection
}
\description{This package uses sophisticated mixed modelling methods
  with the addition of allowing a whole genome approach to detecting
  significant QTL in linkage maps.
}
\details{
\tabular{ll}{
Package: \tab wgaim\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2011-10-24\cr
License: \tab GPL 2\cr
}
Welcome to version 1.x of wgaim! The documentation given in this help file is
only brief and users should eventually consult the soon to be released vignette
that will be a companion to version 1.0 and above.

Version 1.x package highlights:

1) The inclusion of high dimensional genetic data
2) Users have a choice of whole genome marker or interval analysis.
3) Selection and estimation of QTL effects occurs using the
implementation of Verbyla & Taylor (2011).

Deprecated functions:

1) \code{read.interval}: \code{cross2int} is the conversion function
that should be used.
2) \code{wmerge}: This function is not required as the genetic
data is merged within the wgaim call.

This package builds on the \pkg{qtl} package of Broman by including
additional functions for whole genome QTL analysis of a full linkage map
using linear mixed models.

The package provides a user friendly function \code{\link{cross2int}} for the conversion
of "\code{cross}" objects created using \code{read.cross} in Bromans
\pkg{qtl} package into an "\code{interval}" object or use in
\pkg{wgaim}. Specifically, \code{cross2int} performs additional
calculations for imputing missing values on each of the
chromosomes across the full linkage map and also provides users with
genetic distances and recombination fractions for the intervals. The
returned object retains the class structure
of an object created with \code{read.cross} and thus allows further
use with the \pkg{qtl} package if desired.

The package also provides a very neat graphical display of the
chromosomes of a "\code{cross}" object. The method function \code{\link{link.map}} displays
the full or subsetted linkage map according to chromosome or distance as
well as displays non-overlapping marker names on the right hand side.

Modelling of the QTL's is achieved using the functions
\code{\link{wgaim}} which, as its first argument, requires an \code{asreml}
base model. Version 1.0 of \pkg{wgaim} allows users to include high
dimensional genetic components in a \code{wgaim} analysis (See
\code{wgaim.asreml} for more details or the soon to be
released vignette companion to this version). For convenience the
default tracing of results from the asreml models is
outputted to a file for further inspection. For diagnostic purposes, the
outlier statistics from each iteration can be viewed using
\code{\link{out.stat}}. Diagnostics of the likelihood ratio test performed for
each forward step can be displayed using \code{\link{tr.wgaim}}. The function
also displays an incremental probability value matrix of the QTL ascertained at
each forward step of the algorithm.

Summary and print methods are available for the returned "\code{wgaim}"
object and provide users with a detailed report on the QTL, their
size, their flanking markers and significance (including LOD score). The
returned "\code{wgaim}" object may also be plotted using the method
function \code{\link{link.map}}. This function plots the full linkage map
subsetted for chromosome and distance as well as provides shaded
QTL regions and highlighted flanking markers. Plotting of QTL for
multiple traits is also possible (see \code{\link{link.map.default}})

}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis
  Maintainer: Julian Taylor <julian.taylor@csiro.au>
}
\references{
  Verbyla, A. P & Taylor, J. D (2011). High dimensional whole genome
  average interval mapping and a random effects formulation.
  \emph{Theoretical and Applied Genetics}. Submitted.

 Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
 Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
 Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.

 Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
 by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\keyword{ package }
\seealso{
  \code{\link[qtl]{qtl-package}}
}

