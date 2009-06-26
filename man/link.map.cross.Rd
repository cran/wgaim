\name{link.map.cross}
\alias{link.map.cross}
\alias{link.map}
\title{Plot a genetic linkage map}
\description{
Neatly plots the genetic linkage map with marker locations and
marker names. 
}
\usage{
\method{link.map}{cross}(object, chr, max.dist, marker.names = "markers",
tick = FALSE, squash = TRUE, \ldots)
}
\arguments{
  \item{object}{object of class "\code{cross}"}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{max.dist}{a numerical value in cM determining the distance the
    genetic map should be subsetted by}
  \item{marker.names}{a character string naming the type of marker
    information to plot. If "\code{dist}" then distances names
    plotted alongside each chromosome on the left. If "\code{markers}"
    then marker names are plotted instead. Defaults to "\code{markers}"}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names. Defaults to \code{FALSE}}
  \item{squash}{logical value. if \code{TRUE} then creates extra room on
    the left side of the chromosomes. This is useful for plotting trait
    names for QTLs using \code{link.map.wgaim} and \code{link.map.default}}  
  \item{\ldots}{arguments passed to "\code{plot}" to set up the plot
    region. Arguments may also be passed to "\code{text}" for the
    manipulation of the marker names}
}
\details{
This plotting procedure provides a neater visual display of the
chromosomes without marker names overlapping vertically. The plotting region will
adjust itself to ensure that all marker names are in the region. For
this reason the default value for the "\code{cex}" parameter passed to
"\code{text}" should be manipulated until an aesthetic genetic map is reached.

For large maps with many chromosomes, marker names and adjacent
chromosomes will overlap horizontally. For the interest of readability
this has not been corrected. For this particular situation it is
suggested that the user horizontally maximise the plotting window
until no overlapping occurs or subset the genetic map to achieve the
desired result.    
}
\value{This invisibly returns the following list for manipulation with
  \code{link.map.wgaim}
  \item{mt}{A list named by the chromosomes with each element containing the locations of
   the marker names after correcting for overlapping}
  \item{map}{A list named by the chromosomes with each element
    containing the locations of markers on the chromosomes}
  \item{chrpos}{The numerical position of the chromosomes on the
    plotting region} 
}
\author{Julian Taylor}
\seealso{\code{\link{link.map.wgaim}}}
\examples{

data(raccas, package = "wgaim")
link.map(raccas, cex = 0.5) 

}
\keyword{hplot}
