\name{summary.wgaim}
\alias{summary.wgaim}
\alias{print.wgaim}
\alias{print.summary.wgaim}
\title{Summary and print methods for the class "wgaim"}
\description{
Prints a summary of the "\code{wgaim}" object in a presentable format
}
\usage{
\method{summary}{wgaim}(object, parentData, LOD = TRUE, \ldots)
\method{print}{wgaim}(x, parentData, \ldots)
\method{print}{summary.wgaim}(x, \ldots)
}
\arguments{
  \item{object}{an object of class "\code{wgaim}"}
  \item{x}{an object of class "\code{wgaim}" or "\code{summary.wgaim}"
    (see Details)}
  \item{parentData}{A data structure of class "\code{cross}" containing
    the genotypic data as well as the phenotypic data}
  \item{LOD}{Calculate LOD scores for QTL, defaults to \code{TRUE}}
  \item{\ldots}{further arguments passed to or from other methods}
}
\details{
  It is important that the \code{parentData} is not missing in
  \code{summary.wgaim} or \code{print.wgaim} as it
  contains vital summary information about each of the QTL's
  detected. \code{print.summary.wgaim} is just the print method for \code{summary.wgaim}.   
}
\value{A summary of the QTL component of the "\code{wgaim}" object is
    printed to the screen. For each QTL detected, the
    \code{summary.wgaim} prints which chromosome, name and distance of
    each flanking marker, size, z-ratio and LOD
    score. \code{print.wgaim} provides a narrative brief of the QTL's detected.  
}
\references{
  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
  
  Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
  by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 195-211. 
  }
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\seealso{\code{\link{wgaim.asreml}}}
\examples{

\dontrun{
# read in data 

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and merge genotypic with phenotypic

raccasS <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccasM <- wmerge(raccasS, zinc, by = "id")

## base model

zn.fm <- asreml(znconc ~ Type, random = ~ Block + id, data = zinc)

# find QTL's

zn.qtl <- wgaim(zn.fm, parentData = raccasM, trace = "trace.txt",
na.method.X = "include")

# summarise

print(zn.qtl, raccasM)
summary(zn.qtl, raccasM)

}
}
\keyword{regression}

