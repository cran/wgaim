\name{qtlTable}
\alias{qtlTable}
\title{Stack QTL summary information into a super table}
\description{Stack QTL summary information into a super table ready for
  simple exporting
}
\usage{
qtlTable(\ldots, intervalObj = NULL, labels = NULL, columns = "all")
}
\arguments{
  \item{ldots}{list of objects of class "\code{wgaim}". All
    "\code{wgaim}" models must have been analysed with the same genetic
    type (see \code{wgaim.asreml})}
  \item{intervalObj}{a genetic object of class "interval" usually used
    in a \code{wgaim} fit}
  \item{labels}{a vector of character strings determining the trait
    names of each QTL table. if this is NULL then the trait names are
    found from the response of the wgaim model}
  \item{columns}{this can be either a numeric vector determining which columns of the QTL
    summaries should be outputted or "\code{all}" for all
    columns. The default is "\code{all}".}
  \item{\ldots}{a numeric vector determining which columms of the
    summary should be returned (see Details)}
}
\details{
The super table is created by stacking the QTL summaries on top of each
other using the models in \ldots from left to right. An extra column is
created on the left hand side of the stacked table for the trait names
given in the \code{labels} argument. The names are only placed in the
first element of each table with NAs for the rest of the elements. This
then allows simple exporting to spreadsheet packages or with the R/LaTeX
package \pkg{xtable}.
}
\value{A \code{data.frame} object with stacked QTL summaries
 }
 \references{
   Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
   Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
   Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{wgaim}}}
\examples{

\dontrun{

## fit wgaim models

rktgw.qtl <- wgaim(rktgw.asf, phenoData = phenoRxK, intervalObj = genoRxK,
merge.by = "Genotype", trace = "trace.txt", na.method.X = "include")

rkyld.qtl <- wgaim(rkyld.asf, phenoData = phenoRxK, intervalObj = genoRxK,
merge.by = "Genotype", trace = "trace.txt", na.method.X = "include")

## create super table and export

qtlt <- qtlTable(rktgw.qtl, rkyld.qtl, labels = c("Conc.", "Shoot"))
print(xtable(qtlt), file = "superQTL.tex", include.rownames = FALSE)
}
}

\keyword{regression}
