\name{out.stat}
\alias{out.stat}
\title{Plot the interval outlier statistics from
  specififed iteratons of wgaim}
\description{
Plots the interval outlier statistics for
specified iterations of \code{wgaim}. The interval statistics appear as
a trace across the genome separated by chromosome and appropriately
spaced by their distances.
}
\usage{
out.stat(object, intervalObj, iter = NULL, chr =
NULL, \ldots)
}
\arguments{
  \item{object}{object of class "\code{wgaim}"}
  \item{intervalObj}{object of class "\code{cross}" or "\code{interval}"}
  \item{iter}{numeric value determining which iterations will be plotted}
  \item{chr}{character string naming the subset of chromosomes to
    plot. This can only be used when \code{int} is \code{TRUE}}
  \item{\ldots}{arguments passed to "\code{xyplot}" or "\code{barchart}"
  (with some restrictions, see Details)}
}
\details{
By default the interval outlier statistics are plotted in separate
panels for each iteration in a set layout of 5 rows and one column. This
cannot be adjusted and users should not attempt to use the \code{layout}
argument. Viewing multiple pages can be done by specifying the
appropriate iterations using the \code{iter} argument.

The set of QTL are obtained from the model and printed on the plot in their
appropriate positions in each panel.
}
\value{The outlier statistics are plotted in a trellis panel plot.
}
\references{
  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{tr.wgaim}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

zn.qtl <- wgaim(zn.fm, phenoData = zinc, intervalObj = raccasS, trace = "trace.txt",
na.method.X = "include")

# plot QTL interval outlier statistics

out.stat(zn.qtl, raccasM, iter = 1:5, cex = 0.4)

}
}
\keyword{hplot}
