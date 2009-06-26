\name{out.stat}
\alias{out.stat}
\title{Plot the chromosome or interval outlier statistics from
  specififed iteratons of wgaim}
\description{
Plots the outlier statistic for either the chromosome or interval for
specificied iterations of \code{wgaim}. The
chromosome outlier statistics appears as a barchart across the
genome. The interval statistics appear as a trace across the genome
separated by chromosome and appropriately spaced by their distances
}
\usage{
out.stat(object, parentData, int = TRUE, iter = NULL, chr =
NULL, \ldots)
}
\arguments{
  \item{object}{object of class "\code{wgaim}"}
  \item{parentData}{object of class "\code{cross}"}
  \item{int}{logical value. If \code{TRUE} then plot interval
    statistics. If \code{FALSE} then plot chromosome statistics}
  \item{iter}{numeric value determining which iterations will be plotted} 
  \item{chr}{character string naming the subset of chromosomes to
    plot. This can only be used when \code{int} is \code{TRUE}}
  \item{\ldots}{arguments passed to "\code{xyplot}" or "\code{barchart}"
  (with some restrictions, see DETAILS)}
}
\details{
By default the chromosome/interval outlier statistics are plotted in separate
panels for each iteration in a set layout of 5 rows and one column. This
cannot be adjusted and users should not attempt to use the \code{layout}
argument. Viewing multiple pages can be done by specifying the
appropriate iterations using the \code{iter} argument.   

For both the chromosome and interval outlier statistic plots, the set of
QTL are obtained from the model and printed on the plot in their
appropriate positions in each panel.
}
\value{The outlier statistics are plotted in a trellis panel plot. 
}
\author{Julian Taylor}
\seealso{\code{\link{tr.wgaim}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

zn.qtl <- wgaim(zn.fm, parentData = raccasM, trace = "trace.txt",
na.method.X = "include")

# plot QTL interval outlier statistics

out.stat(zn.qtl, raccasM, int = TRUE, iter = 1:5, cex = 0.4)

}
}
\keyword{hplot}
