\name{cross2int}
\alias{cross2int}
\title{Convert a cross genetic object to an interval object}
\description{
Converts an object of class "\code{cross}" to an object with class
"\code{interval}". The function also imputes missing markers.
}
\usage{
cross2int(fullgeno, missgeno = "MartinezCurnow", rem.mark = TRUE,
     id = "id", subset = NULL)
}
\arguments{
  \item{fullgeno}{an object of class "\code{cross}" with restrictions
    (see Details)}
  \item{missgeno}{a character string determining how missing values in
    the linkage map should be imputed. If "\code{Broman}", then missing values are
    imputed according to Bromans rules. If "\code{MartinezCurnow}" then
    missing values are imputed according to the rules of Martinez &
    Curnow (1994) (see reference list). The default is
    "\code{MartinezCurnow}" (see Details).}
  \item{rem.mark}{logical value. If \code{TRUE} redundant markers are deleted
    and placed in the component of the object (see Details). Defaults to
  \code{TRUE}.}
  \item{id}{a character string or name of the unique identifier for each row of genotype
    data (see Details). Defaults to "\code{id}"}
  \item{subset}{a possible character vector naming the subset of
    chromosomes to be returned. Defaults to \code{NULL} implying return
    all chromosomes.}
}
\details{
This function provides the conversion of genetic data objects that have
already been created using \code{read.cross} from Bromans \pkg{qtl}
package to "\code{interval}" objects ready for use with \code{wgaim}.

User should be aware that this function is restricted to populations with
only two genotypes. \code{fullgeno} is checked for the
the class structure \code{c("bc","cross")}. If this is not present an
error is returned.

During the conversion process, missing values are imputed according to
the argument given by \code{missgeno}. This imputation results in a
complete version of the marker data for each chromosome which is then
used to create the interval data "\code{intval}". The complete
marker data for each chromosome can be obtained from the "\code{imputed.data}" element of the
returned list. It is therefore also possible to perform whole genome \emph{marker}
analysis using \code{wgaim}. See \code{wgaim.asreml} for more details.

If \code{rem.mark = TRUE} then markers that are identical are removed before
missing values are imputed. The marker similarity is found by estimating
recombination fractions using \code{est.rf} from the \pkg{qtl} package.
The removed list is returned with the object and placed under
"\code{cor.markers}" for inspection if required.
}
\value{
a list of class "\code{cross}" that also inherits the class
"\code{interval}". The list contains the following components
\item{geno}{
  This is a list with elements named by the corresponding names of the
  chromosomes. Each chromosome is itself a list with six
  elements: "\code{data}" is the actual estimated map matrix with rows
  as individuals named by "\code{id}" and markers as columns; "\code{map}" is a vector of marker positions on the
  corresponding chromosome; "\code{imputed.data}" is identical to "\code{data}"
  matrix but with all NA's replaced by imputed values according to the
  rules of "\code{missgeno}"; "\code{dist}" contains the genetic
  distance between adjacent markers or the genetic distances of the
  intervals; "\code{theta}" contains the recombination fractions for
  each interval; "\code{intval}" contains the recalculated intervals
  based on the recombination fractions and the missing marker information.}
\item{cor.markers}{
If \code{rem.mark = TRUE}, a three column matrix with each row
describing which pairwise markers are correlated and what chromosome
they are from.}
\item{pheno}{
  A data.frame of phenotypic information with rows as individuals read
  in from \code{read.cross}. A copy of the column named by the "\code{id}" argument
  can be found here (see \code{read.cross}).}
}
\references{
  Martinez, O., Curnow. R. N. (1994) Missing markers when estimating
  quantitative trait loci using regression mapping. \emph{Heredity},
  \bold{73}, 198-206.

  Julian Taylor, Arunas Vebyla (2011). R Package wgaim: QTL Analysis in
  Bi-Parental Populations Using Linear Mixed Models. \emph{Journal of
  Statistical Software}, \bold{40}(7), 1-18. URL \url{http://www.jstatsoft.org/v40/i07/}.

  Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
  by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor, Simon Diffey, Ari Verblya and Brian Cullis}
\seealso{\code{\link[qtl]{read.cross}}}
\examples{

\dontrun{
# read in linkage map from a rotated .CSV file with "id" as the
# identifier for each unique row

wgpath <- system.file("extdata", package = "wgaim")
racca <- read.cross("csvr", file="raccas.csv", genotypes=c("AA","AB"),
         na.strings = c("-", "NA"), dir = wgpath)
raccas <- cross2int(racca, missgeno="MartinezCurnow", id = "id")

# plot linkage map

link.map(raccas, cex = 0.5)

}
}
\keyword{regression}