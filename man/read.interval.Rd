\name{read.interval}
\alias{read.interval}
\title{Read in and estimate a linkage map and interval information}
\description{
Reads genotypic marker data in a wide range of formats and estimates the
appropriate linkage map via Bromans \code{read.cross} function (see
Details). Additionally, this function imputes missing markers and also
estimates interval information appropriate for further analysis.
}
\usage{
read.interval(format = c("csv", "csvr", "csvs", "csvsr", "mm",
"qtx", "tlcart", "gary", "karl"), dir = getwd(), file, genfile, mapfile,
phefile, chridfile, mnamesfile, pnamesfile, na.strings = c("-", "NA"),
genotypes = c("A", "B"), estimate.map = TRUE, convertXdata = TRUE,
missgeno = "MartinezCurnow", rem.mark = TRUE, id = "id",
subset = NULL, \ldots)
}
\arguments{
  \item{format}{see \code{read.cross}}
  \item{dir}{see \code{read.cross}}
  \item{file}{see \code{read.cross}}
  \item{genfile}{see \code{read.cross}}
  \item{mapfile}{see \code{read.cross}}
  \item{phefile}{see \code{read.cross}}
  \item{chridfile}{see \code{read.cross}}
  \item{mnamesfile}{see \code{read.cross}}
  \item{pnamesfile}{see \code{read.cross}}
  \item{na.strings}{see \code{read.cross}}
  \item{genotypes}{see \code{read.cross}}
  \item{estimate.map}{see \code{read.cross}}
  \item{convertXdata}{see \code{read.cross}} 
  \item{missgeno}{a character string determining how missing values in
    the linkage map should be imputed. If "\code{Broman}", then missing values are
    imputed according to Bromans rules. If "\code{MartinezCurnow}" then
    missing values are imputed according to the rules of Martinez &
    Curnow (1994) (see reference list). The default is "\code{MartinezCurnow}"}
  \item{rem.mark}{logical value. If \code{TRUE} redundant markers are deleted
    and placed in the component of the object (see details). Defaults to
  \code{TRUE}.}
  \item{id}{the name of the unique identifier for each row of genotype
    data (see details). Defaults to "\code{id}"}
  \item{subset}{A possible character vector naming the subset of
    chromosomes to be returned. Defaults to \code{NULL} implying return
    all chromosomes.}
  \item{\ldots}{Any other arguments to be passed to \code{read.cross} for
  the appropriate reading in of the data from a file (see
  \code{read.cross} for more details).}
}
\details{
This function uses a *modified* verison of Bromans' \code{read.cross}
function to read in the genotypic data and estimate the linkage
map. This *strict* modification only allows two genotypes which are the
defaults "A" and "B". Attempting more than this will cause
an error.

Bromans' \code{read.cross} and this function by default read genotypic information as
well as phenotypic information if the data contains any. In general,
this phenotypic component of the file *should* contain, at least, a column
to identify each row of the data uniquely. The argument \code{id} is
used for the purpose of identifying this column. An error will be
induced otherwise.

The output of this function differs slightly from \code{read.cross} by
the inclusion of additional components to the returned object. To ensure
that functions of the \pkg{qtl} package could still be used with the returned
object from this function the class structure from \code{read.cross} has
been maintained. The returned object also inherits the class
"\code{interval}" for forwards compatibility with functions available in
the \code{wgaim} package. 
}
\value{
a list of class "\code{cross}" that also inherits the class
"\code{interval}". The list contains the following components
\item{geno}{
  This is a list with elements named by the corresponding names of the
  chromosomes. Each chromosome is itself a list with six
  elements: "\code{data}" is the actual estimated map matrix with rows
  as individuals named by "\code{id}" and markers as columns; "\code{map}" is a vector of marker positions on the
  corresponding chromosome; "\code{argmax}" is identical to "\code{data}"
  matrix but with all NA's replaced by imputed values according to the
  rules of "\code{missgeno}"; "\code{dist}" contains the genetic
  distance between adjacent markers or the genetic distances of the
  intervals; "\code{theta}" contains the recombination fractions for
  each interval; "\code{interval}" contains the recalculated intervals
  based on the recombination fractions and the missing marker information.}
\item{cor.markers}{
If \code{rem.mark = TRUE}, a three column matrix with each row
describing which pairwise markers are correlated and what chromosome
they are from.}
\item{pheno}{
  A data.frame of phenotypic information with rows as individuals read
  in from \code{read.cross}. A copy of the column named by the "\code{id}" argument
  can be found here (see \code{read.cross})}
}  
\references{
  Martinez, O., Curnow. R. N. (1994) Missing markers when estimating
  quantitative trait loci using regression mapping. \emph{Heredity},
  \bold{73}, 198-206.

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

racca <- read.interval("csvr", file="raccasgroups.csv", genotypes=c("A","B"),
missgeno="MartinezCurnow", id = "id", na.strings = c("-", "NA"))

# plot linkage map

link.map(racca, cex = 0.5)
}
}
\keyword{regression}
