\name{link.map.default}
\alias{link.map.default}
\title{Plot a genetic linkage map with QTL's for multiple traits}
\description{
Neatly plots the genetic linkage map with marker locations,
marker names and highlights QTL's with their associated flanking markers
for multiple traits obtained from a list of \code{wgaim} fits.
}
\usage{
\method{link.map}{default}(object, parentData, chr, max.dist, marker.names = TRUE, 
    clist = list(qcol = rainbow(length(object)), mcol = "red", 
        tcol = rainbow(length(object))), trait.labels = NULL, tick = FALSE, \ldots)
}
\arguments{
  \item{object}{a list object with elements inheriting the class "\code{wgaim}"}
  \item{parentData}{object of class "\code{cross}"}
  \item{chr}{character string naming the subset of chromosomes to plot}
  \item{max.dist}{a numerical value in cM determining the distance the
    genetic map should be subsetted by}
  \item{marker.names}{logical value. If \code{TRUE} then marker names
    are plotted alongside each chromosome on the left. Defaults to
    \code{TRUE}}
  \item{clist}{named list of colors used to highlight the QTL regions and
    their flanking markers. \code{qcol} is the colors of the QTL
    regions (defaults to \code{rainbow(n)} where \code{n} is the length
    of \code{object}). \code{mcol} is the color the flanking markers. \code{tcol}
    is the color of the trait names used in each model (defaults to the
    same color as the QTL regions). See \code{par} for color options} 
  \item{trait.labels}{character string naming the trait used in the
    model object, defaults to the names of the traits used in each
    model.}
  \item{tick}{logical value. If \code{TRUE} then an axis with tick marks are
    generated for the chromosome names}
  \item{\ldots}{arguments passed to "\code{plot}" to set up the plot
    region. Arguments may also be passed to "\code{text}" for the
    manipulation of the marker names}
}
\details{
This plotting procedure is a wrapper for \code{link.map.wgaim} and displays
QTLs for multiple traits obtained from a list of models given by \code{object}. 
Alternative labels for the traits can be given, in model order, using
\code{trait.labels}.

Color specific highligting of the QTLs is also available using
\code{clist}. This differs slightly from \code{link.map.wgaim}. Here
the \code{qcol} and \code{tcol} should be given a set of colors equal to
the length of \code{object}. Let \code{n} be the length of \code{object}.
Then if \code{qcol} is \code{NULL} or length of \code{qcol} is not equal
to \code{n} then it defaults to \code{rainbow(n)}. If \cite{tcol} is
\code{NULL} or length of \code{tcol} is not equal to \code{n} or 1 then
it defaults to the colors of \code{qcol}. Examples of different color
combinations are given below.
}
\value{The genetic linkage map is plotted with shaded QTL regions and
  highlighted flanking markers.
}
\author{Julian Taylor}
\seealso{\code{\link{link.map.cross}}, \code{\link{link.map.wgaim}}}
\examples{

\dontrun{
## fit wgaim models

zn.qtl <- wgaim(zn.fm, parentData = raccasM, trace = "trace.txt",
na.method.X = "include")

zn.qtlS <- wgaim(zn.fmS, parentData = raccasM, trace = "trace.txt",
na.method.X = "include")

## plot QTL intervals

# matching rainbow QTL color and trait names, red flanking markers
# (default) and gray background markers.

link.map(list(zn.qtl,zn.qtlS), raccasM, cex = 0.6, col = "gray")

# rainbow QTL color and black trait names, red flanking markers
# (default) and gray background markers.

link.map(list(zn.qtl,zn.qtlS), raccasM, cex = 0.6, clist(tcol = "black",
mcol = "red"), col = "gray")

# monochromatic plot: gray QTLs, black trait names, black flanking
# markers and gray background markers

link.map(list(zn.qtl,zn.qtlS), raccasM, cex = 0.6, clist(qcol =
rep(gray(0.8), 2), tcol = "black", mcol = "black"), col = "gray")

}
}
\keyword{hplot}