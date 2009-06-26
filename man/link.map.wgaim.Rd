\name{link.map.wgaim}
\alias{link.map.wgaim}
\title{Plot a genetic linkage map with QTL's}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL's with their associated flanking markers
obtained from a fit to \code{wgaim}.
}
\usage{
\method{link.map}{wgaim}(object, parentData, chr, max.dist, marker.names
   = "markers", clist = list(qcol = "light blue", mcol = "red",
   tcol = "light blue"), trait.labels = NULL, tick = FALSE, \ldots)
}
\arguments{
  \item{object}{object of class "\code{wgaim}"}
  \item{parentData}{object of class "\code{cross}"}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{max.dist}{a numerical value in cM determining the distance the
    genetic map should be subsetted by}
   \item{marker.names}{a character string naming the type of marker
    information to plot. If "\code{dist}" then distances names
    plotted alongside each chromosome on the left. If "\code{markers}"
    then marker names are plotted instead. Defaults to "\code{markers}".}
  \item{clist}{named list of colours used to highlight the QTL regions and
    their flanking markers. \code{qcol} is the color of the QTL
    regions. \code{mcol} is the color the flanking markers. \code{tcol}
    is the color of the trait name used in the model object (see \code{par} for colour options)} 
  \item{trait.labels}{character string naming the trait used in the
    model object}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names}
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

The \code{codeList} arguments \code{qcol}, \code{mcol} and \code{tcol} have been added for personal
colour highlighting of the QTL regions, flanking markers and trait names. For greater
flexibility the procedure may also be given the usual \code{col}
argument that will be passed to the other markers.

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

link.map(zn.qtl, raccasM, cex = 0.6, clist = list(mcol = "red"),
     col = "gray")

}
}
\keyword{hplot}
