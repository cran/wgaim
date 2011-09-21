\name{link.map.wgaim}
\alias{link.map.wgaim}
\title{Plot a genetic linkage map with QTL's}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL's with their associated flanking markers
obtained from a fit to \code{wgaim}.
}
\usage{
\method{link.map}{wgaim}(object, intervalObj, chr, max.dist,
    marker.names = "markers", list.col = list(q.col = "light blue",
    m.col = "red", t.col = "light blue"), list.cex = list(t.cex = 0.6,
    m.cex = 0.6), trait.labels = NULL, tick = FALSE, \ldots)
}
\arguments{
  \item{object}{object of class "\code{wgaim}"}
  \item{intervalObj}{object of class "\code{cross}" or "\code{interval}"}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{max.dist}{a numerical value in cM determining the distance the
    genetic map should be subsetted by}
   \item{marker.names}{a character string naming the type of marker
    information to plot. If "\code{dist}" then distances names
    plotted alongside each chromosome on the left. If "\code{markers}"
    then marker names are plotted instead. Defaults to "\code{markers}".}
  \item{list.col}{named list of colours used to highlight the QTL regions and
    their flanking markers. \code{q.col} is the color of the QTL
    regions. \code{m.col} is the color the flanking markers. \code{t.col}
    is the color of the trait name used in the model object (see
    \code{par} for colour options)}
  \item{list.cex}{a named list object containing the character expansion
    factors for the marker names \code{m.cex} and the trait labels \code{t.cex}}
  \item{trait.labels}{character string naming the trait used in the
    model object}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names}
  \item{\ldots}{arguments passed to "\code{plot}" to set up the plot
    region and plot any sybols if required}
}
\details{
This plotting procedure builds on \code{link.map.cross} by adding the
QTL regions to the map and highlighting the appropriate markers obtained
from a fit to \code{wgaim}. If the linkage map is subsetted and QTL
regions fall outside the remaining map a warning will be given that
the QTL have been omitted from the display.

The \code{list.col} arguments \code{q.col}, \code{m.col} and \code{t.col} have been added for personal
colour highlighting of the QTL regions, flanking markers and trait names. For greater
flexibility the procedure may also be given the usual \code{col}
argument that will be passed to the other markers.

The \code{list.cex} argument can be used to manipulate the character expansion of
the marker names using \code{m.cex} or the character expansion of the
\code{trait.labels} using \code{t.cex}. If a "marker" analysis has been
performed then \code{pch} is used to plot a symbol at the
location of the QTL. This character can be changed using the usual
arguments such as \code{pch} or \code{cex} that are passed through the
\ldots argument.
}
\value{For an "interval" analysis, the genetic linkage map is plotted with shaded QTL regions and
  highlighted flanking markers. For a "marker" analysis, a symbol is
  placed at the QTL locations and the markers are highlighted.
}
\references{
  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{link.map.cross}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

zn.qtl <- wgaim(zn.fm, phenoData = zinc, intervalObj = raccas,
merge.by = "id", trace = "trace.txt", na.method.X = "include")

# plot QTL

link.map(zn.qtl, raccas, list.col = list(m.col = "red"), col = "gray")

}
}
\keyword{hplot}
