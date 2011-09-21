\name{wgaim.asreml}
\alias{wgaim.asreml}
\alias{wgaim}
\title{wgaim method for class "asreml"}
\description{Fits an iterative Whole Genome Average Interval Mapping
  (wgaim) model for QTL detection
}
\usage{
\method{wgaim}{asreml}(baseModel, phenoData, intervalObj, merge.by = NULL,
         gen.type = "interval", method = "random", TypeI = 0.05, attempts = 5,
         trace = TRUE, verboseLev = 0, \ldots)
}
\arguments{
  \item{baseModel}{a model object of class "\code{asreml}" usually representing
    a base model with which to build the qtl model.
}
  \item{phenoData}{a data frame containing the phenotypic
    elements used to fit \code{baseModel}. This data is checked against the base
    models data
}
  \item{intervalObj}{a list object containing the genotypic data, usually an "\code{interval}"
    object obtained from using \code{cross2int}. This object may contain
    many more markers than observations (see Details).
}
 \item{merge.by}{a character string or name of the column(s) in \code{phenoData} and
   \code{intervalObj} to merge the phenotypic and genotypic data sets.
 }
\item{gen.type}{a character string determining the type of genetic data to
  be used in the analysis. Possibilities are "\code{interval}" and
  "\code{markers}". The default is "\code{interval}". (see Details).
}
\item{method}{a character string determining the type of algortithm to
  be used in the analysis. Possibilities are "\code{random}" and
  "\code{fixed}". The default is "\code{random}". (see Details).
}
  \item{TypeI}{a numerical value determining the level of significance
    for detecting a QTL. The default is 0.05.
  }
  \item{attempts}{An integer representing the number of attempts at convergence for the fixed or
    random qtl model. The default is 5.}
  \item{trace}{An automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace machanisms will display a
  message if a QTL is detected.}
 \item{verboseLev}{numerical value, either 0 or 1, determining the level of tracing
   outputted during execution of the algorithm A 0 value will produce
   the standard model fitting output from the fitted ASReml models
   involved in the forward selection. A value of 1 will add a table of
   Chromosome and Interval outlier statistics for each iteration}
  \item{\ldots}{Any other extra arguments to be passed to each of the
    \code{asreml} calls. These may also include \code{asreml.control} arguments.}
}
\details{
In the initial call to \code{wgaim.asreml}, the marker or interval
information is collated from \code{intervalObj}. If \code{gen.type =
  "interval"} then midpoints of intervals are collated from the
"\code{intval}" component of \code{intervalObj}. If \code{gen.type =
  "markers"} then markers are collated from the "\code{imputed.data}"
component of \code{intervalObj} (It should be noted that a
"\code{marker}" analysis is less efficient than an "\code{interval}" analysis as
it does not take into account the correlation of the marker effects in
the specificity of the model; see Verbyla et. al, 2007).

The \code{method} argument in \code{wgaim.asreml} allows the user access
to two algorithms. If \code{method = "fixed"} then the forward
selecction algorithm uses the legacy code packaged with pre-1.0 versions
of \pkg{wgaim}. This code uses the chromosome statistics as a guide to
its selection of QTL. This selection process is now known to be flawed
when small linkage groups are present. For this reason it is suggested
that this option only be used when there are a moderate number of
markers on each linkage group.

This version of \pkg{wgaim} allows high dimensional marker information
to be analysed. A simple transformation of the collated high dimensional
marker set shows that it may be reduced to the number of genetic lines
used in the analysis. This transformation is internal to the
\code{wgaim.asreml} call and users can now expect a considerably large
acceleration in the performance of \pkg{wgaim}.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.
}
\value{An object of class "\code{wgaim}" which also inherits the class
  "\code{asreml}" by default. The object returned is actually an \code{asreml}
  object (see \code{asreml.object}) with the addition of components from
  the QTL detection listed below.

  \item{QTL}{A list of components from the significant QTL detected
    including a character vector of the significant QTL along with a
    vector of the QTL effect sizes. There are also a number of diagnostic meausres that
    can be found in \code{diag} that are used in conjunction with
    \code{tr.wgaim} and \code{out.stat}.}
 }
 \references{
   Verbyla, A. P & Taylor, J. D (2011). High dimensional whole genome
   average interval mapping and a random effects formulation.
   \emph{Theoretical and Applied Genetics}. Submitted.

   Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
   Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
   Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.

 Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
 by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\seealso{\code{\link{print.wgaim}}, \code{\link{summary.wgaim}}}
\examples{

\dontrun{
# read in data

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and convert to "interval" object

raccas <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccas <- cross2int(raccas, missgeno = "Martinez")

# base model

zn.fm <- asreml(znconc ~ Type, random = ~ Block + id, data = zinc)

# find QTL's

zn.qtl <- wgaim(zn.fm, phenoData = zinc, intervalObj = raccas,
merge.by = "id", trace = "trace.txt", na.method.X = "include")
}

}
\keyword{regression}
