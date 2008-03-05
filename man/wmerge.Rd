\name{wmerge}
\alias{wmerge}
\title{Merge phenotypic data with genotypic data}
\description{ Merge phenotypic data with genotypic data from an object
  of class "\code{interval}" 
}
\usage{
wmerge(geno, pheno, by = NULL, \ldots)
}
\arguments{
  \item{geno}{genotypic data from an object of class
    "\code{interval}" (see \code{read.interval}}
  \item{pheno}{phenotypic data which may be give as a \code{data.frame}
    or a file}
  \item{by}{character string defining the name of the column(s) with
    which to merge the phenotypic data with the genotypic data}
  \item{\ldots}{arguments passed to \code{read.table} if \code{pheno} is
    a filename}
}
\details{
The \code{\ldots} argument actually passes extra arguments to
\code{asreml.read.table} which in turn passes identical arguments to
\code{read.table}. Therefore if \code{header} or \code{col.names} is
set, names of columns with a capital letter are converted to a factor
(see \code{asreml.read.table}.)   
  
This function provides a fail safe mechanism with which to merge
large scale genotypic data with phenotypic data. It is important that
both data sets contain a matching column named with the argument
"\code{by}".

As the phenotypic data is inherently larger in size than the genotypic
data the merging process has been simplified to ensure that the correct
combined data set is returned. Any names of the genotypic "\code{by}"
column that do not appear in the names of the phenotypic "\code{by}" column
are dropped from the genotypic data as there is no phenotypic
information for these uniquely named rows. Any names of the phenotypic "\code{by}"
column that do not appear in the names of the genotypic "\code{by}"
column induce a row of NA's to be placed in the genotype data. This is
to ensure that the phenotypic data (which may include columns containing
important design information) remains intact.     
}
\value{a list contaning identical components as returned by
  \code{wgaim.read.cross} or possibly \code{read.cross} plus the
  additional following components
  \item{pheno.dat}{A copy of the phenotypic data used in the merge}
  \item{full.data}{A copy of the full data, i.e the genotypic data merged
    with the phenotypic data using the matching column named by the
    argument "\code{by}"} 
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\seealso{\code{\link{wgaim}}, \code{\link{read.interval}}}
\examples{

\dontrun{

# read in data 

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and merge genotypic with phenotypic

raccasS <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccasM <- wmerge(raccasS, zinc, by = "id")

}
}

\keyword{regression}

