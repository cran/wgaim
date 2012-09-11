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
  Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
  dimensional random whole genome average (QTL) interval mapping approach.
  \emph{Genetics Research}. Submitted.

   Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
   Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
   Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.
}
\author{Julian Taylor}
\seealso{\code{\link{wgaim}}}
\examples{

\dontrun{
# read in data

data(phenoRxK, package = "wgaim")
data(genoRxK, package = "wgaim")

# subset linkage map and convert to "interval" object

genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))
genoRxK <- cross2int(genoRxK, missgeno = "Martinez",
                    id = "Genotype", map.function = "kosambi")

# base model

yield.fm <- asreml(yield ~ Type, random = ~ Genotype + Row + Range,
                   rcov = ~ ar1(Range):ar1(Row), data = phenoRxK)

# find QTL

yield.qtl <- wgaim(yield.fm, phenoData = phenoRxK, intervalObj = genoRxK,
   merge.by = "Genotype", gen.type = "interval", method = "fixed",
   selection = "interval", trace = "trace.txt", na.method.X = "include")

# diagnostic check

tr(yield.qtl, digits = 4)
}
}

\keyword{regression}

