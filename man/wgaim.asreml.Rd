\name{wgaim.asreml}
\alias{wgaim.asreml}
\alias{wgaim}
\title{wgaim method for class "asreml"}
\description{Fits an iterative Whole Genome Average Interval Mapping
  (wgaim) model for QTL detection
}
\usage{
\method{wgaim}{asreml}(baseModel, phenoData, intervalObj, merge.by = NULL,
         gen.type = "interval", method = "fixed", selection = "interval",
         breakout = -1, TypeI = 0.05, attempts = 5, trace = TRUE,
         verboseLev = 0, \ldots)
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
\item{method}{a character string determining the type of algorithm to
  be used in the analysis. Possibilities are "\code{random}" and
  "\code{fixed}". The default is "\code{random}". (see Details).
}
\item{selection}{a character string determining the type of selection
  method that is used to select QTL in the analysis. Possibilities are "\code{interval}" and
  "\code{chromosome}". The default is "\code{interval}". (see Details).
}
\item{breakout}{A numerical integer equivalent to the number of
  iterations the algorithm performs before breaking out. The default is
  -1 which ensures the algorithm finds all QTL before halting. (see
  Details)
}
\item{TypeI}{a numerical value determining the level of significance
    for detecting a QTL. The default is 0.05.
  }
\item{attempts}{An integer representing the number of attempts at convergence for the fixed or
  random QTL model. The default is 5.
}
\item{trace}{An automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace machanisms will display a
  message if a QTL is detected.
}
 \item{verboseLev}{numerical value, either 0 or 1, determining the level of tracing
   outputted during execution of the algorithm A 0 value will produce
   the standard model fitting output from the fitted ASReml models
   involved in the forward selection. A value of 1 will add a table of
   chromosome and interval outlier statistics for each iteration.
}
\item{\ldots}{Any other extra arguments to be passed to each of the
    \code{asreml} calls. These may also include \code{asreml.control}
    arguments.
}
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
the specificity of the model; see Verbyla et. al (2007).

The \code{method} argument in \code{wgaim.asreml} allows the user access
to two algorithms. If \code{method = "fixed"} the algorithm
places selected QTL as an additive set of fixed
effects in the model as the forward selection algorithm proceeds. If \code{method = "random"}
places selected QTL in the random part of the model as
an additive set of random effects. This new formulation is outlined in
Verbyla et. al (2012).

The \code{selection} argument determines the type of selection algorithm
for the analysis. If \code{selection = "chromosome"} then outlier
statistics for each chromosome are calculated and the largest
chromosome or linkage group is chosen. The largest marker/interval outlier
statistic in this linkage group is then selected as the putative QTL. If
\code{selection = "interval"}, only marker/interval statistics are calculated
and the largest marker/interval is chosen as the putative QTL.

Note: If a genetic map has a small number of markers on a linkage group
then using \code{selection = "chromosome"} as the selection algorithm
is known to be flawed (see Verbyla & Taylor, 2012). For this reason it is suggested
that this option only be used when there are a moderate number of
markers on each linkage group.

Users can now break out of the algorithm using the \code{breakout}
argument. If a numerical value greater than zero then the forward
selection algorithm breaks at the iteration equal to that value and
returns the collected QTL information to this point. This includes
fixed/random QTL effects, diaganostic components such as interval/marker
BLUPs and outlier statistics as well as the trace components of the
algorithm. It should be noted that the algorithm breaks out before QTL
have been moved to the fixed/random effects and estimated. Therefore a
value of 1 will not return any estimated QTL effects but it will return
diagnostic information of the first whole genome random effect model.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.

This version of \pkg{wgaim} allows high dimensional marker information
to be analysed. A simple transformation of the collated high dimensional
marker set shows that it may be reduced to the number of genetic lines
used in the analysis. This transformation is internal to the
\code{wgaim.asreml} call and users can now expect a considerably large
acceleration in the performance of \pkg{wgaim}.
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
   Verbyla, A. P & Taylor, J. D, Verbyla, K. L (2012). RWGAIM: An efficient high
  dimensional random whole genome average (QTL) interval mapping approach.
  \emph{Genetics Research}. \bold{94}, 291-306.

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

data(phenoRxK, package = "wgaim")
data(genoRxK, package = "wgaim")

# subset linkage map and convert to "interval" object

genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))
genoRxK <- cross2int(genoRxK, missgeno = "Martinez",
                    id = "Genotype", map.function = "kosambi")

# base model

rkyld.asf <- asreml(yld ~ Type + lrow, random = ~ Genotype + Range,
                   rcov = ~ ar1(Range):ar1(Row), data = phenoRxK)

# find QTL

rkyld.qtl <- wgaim(rkyld.asf, phenoData = phenoRxK, intervalObj = genoRxK,
   merge.by = "Genotype", gen.type = "interval", method = "fixed",
   selection = "interval", trace = "trace.txt", na.method.X = "include")

}
}
\keyword{regression}
