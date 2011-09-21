\name{tr.wgaim}
\alias{tr.wgaim}
\title{Display diagnostic information about wgaim QTL model}
\description{Displays diagnostic infomation about QTL detection and
  significance for the sequence of models used in a \code{wgaim} fit.
}
\usage{
\method{tr}{wgaim}(object, iter = 1:length(object$QTL$effects),
      diag.out = TRUE, \ldots)
}
\arguments{
  \item{object}{an object of class "\code{wgaim}"}
  \item{iter}{a vector of integers describing the rows of the p-value matrix to display}
  \item{diag.out}{logical value. If \code{TRUE} then diagnostic information
about the testing of the genetic variance is given for all iterations.}
  \item{\ldots}{arguments passed to \code{print.default} for displaying
of information}
}
\details{
By default the printing of the objects occur with arguments \code{quote = FALSE}
and \code{right = TRUE}. Users should avoid including these arguments.
}
\value{A probability value matrix for the successvie QTL selected is displayed
with rows according to the iterations specified. If \code{diag.out =
TRUE} then a matrix with rows consisting of the likelihood with
genetic variance, the likelihood without genetic variance (NULL model),
the statistic and the p-value for the statistic.
 }
 \references{
   Verbyla, A. P & Taylor, J. D (2011). High dimensional whole genome
   average interval mapping and a random effects formulation.
   \emph{Theoretical and Applied Genetics}. Submitted.

   Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
   Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
   Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{wgaim}}}
\examples{

\dontrun{
# read in data

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and convert to "interval" object

raccasS <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccasS <- cross2int(raccasS, missgeno = "Martinez")

# base model

zn.fm <- asreml(znconc ~ Type, random = ~ Block + id, data = zinc)

# find QTL's

zn.qtl <- wgaim(zn.fm, phenoData = zinc, intervalObj = raccasS, trace = "trace.txt",
na.method.X = "include")

# diagnostic check

tr(zn.qtl, digits = 4)
}
}

\keyword{regression}

