\name{link.map.wgaim}
\alias{link.map.wgaim}
\title{Plot a genetic linkage map with QTL's}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL's with their associated flanking markers
obtained from a fit to \code{wgaim}.
}
\usage{
\method{link.map}{wgaim}(object, parentData, chr, max.dist, marker.names = TRUE,
qcol = "light blue", mcol = "red", \ldots)
}
\arguments{
  \item{object}{object of class "\code{wgaim}"}
  \item{parentData}{object of class "\code{cross}"}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{max.dist}{a numerical value in cM determining the distance the
    genetic map should be subsetted by}
  \item{marker.names}{logical value. If \code{TRUE} then marker names
    are plotted alongside each chromosome on the left. Defaults to \code{TRUE}}
  \item{qcol}{color of QTL regions (see \code{par} for colour options)}
  \item{mcol}{color of QTL flanking markers (see \code{par} for colour options)}
  \item{\ldots}{arguments passed to "\code{plot}" to set up the plot
    region. Arguments may also be passed to "\code{text}" for the
    manipulation of the marker names}
}
\details{
This plotting procedure builds on \code{link.map.cross} by adding the
QTL regions to the map and highlighting the flanking markers obtained
from a fit to \code{wgaim}. If the linkage map is subsetted and QTL
regions fall outside the remaining map a warning will be given that
the QTL's have been omitted from the display.

The arguments \code{mcol} and \code{qcol} have been added for personal
colour highlighting of the QTL regions and flanking markers. For greater
flexibility the procedure may also be given the usual \code{col}
argument that will be passed to the other markers.

At the current stage the plotting of the linkage map only occurs
horizontally.
}
\value{The genetic linkage map is plotted with shaded QTL regions and
  highlighted flanking markers.
}
\author{Julian Taylor}
\seealso{\code{\link{link.map.cross}}, \code{\link{wgaim}}}
\examples{

\dontrun{
# fit wgaim model

zn.qtl <- wgaim(zn.fm, parentData = raccasM, trace = "trace.txt",
na.method.X = "include")

# plot QTL intervals

link.map(zn.qtl, raccasM, cex = 0.6, mcol = "red", col = "blue")

}
}
\keyword{hplot}
