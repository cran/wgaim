\name{summary.wgaim}
\alias{summary.wgaim}
\alias{print.wgaim}
\alias{print.summary.wgaim}
\title{Summary and print methods for the class "wgaim"}
\description{
Prints a summary of the "\code{wgaim}" object in a presentable format
}
\usage{
\method{summary}{wgaim}(object, intervalObj, LOD = TRUE, \ldots)
\method{print}{wgaim}(x, intervalObj, \ldots)
}
\arguments{
  \item{object}{an object of class "\code{wgaim}"}
  \item{x}{an object of class "\code{wgaim}" (see Details)}
  \item{intervalObj}{a data structure of class "\code{cross}" or
    "\code{interval}" containing the genotypic data}
  \item{LOD}{logical value. If TRUE LOD scores for QTL are calculated, defaults to \code{TRUE}}
  \item{\ldots}{further arguments passed to or from other methods}
}
\details{
  It is important that the \code{intervalObj} is not missing in
  \code{summary.wgaim} or \code{print.wgaim} as it
  contains vital summary information about each of the QTL
  detected.

  As the WGAIM algorithm now places the selected QTL effects in a
  separate random component term of the model, they are summarised
  appropriately using a probablistic argument based on the conditional
  distribution of the QTL sizes given the data (see Verbyla & Taylor,
  2011 in References) Thus, for each QTL, a value is
  calculated that represents the probability that the QTL size is greater
  than zero (or less than zero if the effect is negative).
}
\value{A summary of the QTL component of the "\code{wgaim}" object is
    printed to the screen. For each QTL detected, if an "interval"
    analysis was performed then \code{summary.wgaim} prints which
    chromosome, name and distance of each flanking marker, size, probability and LOD
    score (based on the random coefficient) if desired. If a "marker" analysis was performed then the
    chromosome, name and distance of the associated marker, size, probability
    and LOD score are printed. \code{print.wgaim} provides a narrative brief of the QTL's detected.
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
    Applied Genetics}, \bold{116}, 195-211.
  }
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\seealso{\code{\link{wgaim.asreml}}}
\examples{

\dontrun{
# read in data

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and convert to "interval" object

raccasS <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccasS <- cross2int(raccasS, missgeno = "Martinez")

## base model

zn.fm <- asreml(znconc ~ Type, random = ~ Block + id, data = zinc)

# find QTL

zn.qtl <- wgaim(zn.fm, phenoData = zinc, intervalObj = raccasS, trace = "trace.txt",
na.method.X = "include")

# summarise

print(zn.qtl, raccasM)
summary(zn.qtl, raccasM)

}
}
\keyword{regression}

