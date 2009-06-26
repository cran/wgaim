\name{wgaim-package}
\alias{wgaim-package}
\docType{package}
\title{Whole Genome Average Interval Mapping (wgaim) for QTL detection
}
\description{This package uses sophisticated mixed modelling methods
  with the addition of allowing a whole genome approach to detecting
  significant QTL's in linkage maps.
}
\details{  
\tabular{ll}{
Package: \tab wgaim\cr
Type: \tab Package\cr
Version: \tab 0.99-1\cr
Date: \tab 2008-03-05\cr
License: \tab GPL 2\cr
}
This package builds on the \pkg{qtl} package of Broman by including
additional functions for the modelling and summary of QTL intervals from
the full linkage map.

The package provides user friendly functions for the reading in and
combining of genotypic and phenorypic data. Specifically, users can read
read in genotypic data using \code{\link{read.interval}} which is a modified
version of Bromans \code{read.cross}. \code{\link{read.interval}}
performs additional calculations for imputing missing values on each of the
chromosomes across the full linkage map and also provides users with
genetic distances and recombination fractions for the intervals. There
is also a function \code{\link{cross2int}} that provides the conversion
of "\code{cross}" objects created from the \pkg{qtl} to
"\code{interval}" objects for use in \pkg{wgaim}. Using
\code{\link{wmerge}} the phenotypic data can then be easily
merged with the estimated map produced from \code{\link{read.interval}}
or \code{\link{cross2int}}. The returned object retains the class structure
of an object created with \code{read.cross} and thus allows further
use with the \pkg{qtl} package if desired. 

The package also provides a very neat graphical display of the
chromosomes of a "\code{cross}" object. The method function \code{\link{link.map}} displays
the full or subsetted linkage map according to chromosome or distance as
well as displays non-overlapping marker names on the left hand side.

Modelling of the QTL's is achieved using the functions
\code{\link{wgaim}} which, as its first argument, requires an \code{asreml}
base model. As the WGAIM method requires iterative fitting of
\code{asreml} models users should be patient for results. For
convenience the default tracing of results from the asreml models is
outputted to a file for further inspection. For diagnostic purposes, the
outlier statistics from each iteration can be viewed using
\code{out.stat}. Diagnostics of the likelihood ratio test performed for
each forward step can be displayed using \code{tr.wgaim}. The function
also displays an incremental p-value matrix of the QTL ascertained at
each forward step of the algorithm.   

Summary and print methods are available for the returned "\code{wgaim}"
object and provide users with a detailed report on the QTL, their
size, their flanking markers and significance (including LOD score). The
returned "\code{wgaim}" object may also be plotted using the method
function \code{\link{link.map}}. This function plots the full linkage map
subsetted for chromosome and distance as well as provides shaded
QTL regions and highlighted flanking markers.
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis
  Maintainer: Julian Taylor <julian.taylor@csiro.au>
}
\references{
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

