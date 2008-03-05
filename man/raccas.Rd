\name{raccas}
\alias{raccas}
\docType{data}
\title{Genotypic marker data}
\description{Genotypic marker data represeting a full linkage map from a
  Double Haploid population
}
\usage{data(raccas)}
\format{A data object of class "\code{cross}" inheriting the class
  "\code{interval}".
}
\details{This is data on 93 individuals from a Double Haploid population
  typed at 663 markers over 40 chromosomes. Coincident markers have been
  omitted and the data set has been reduced to 468 markers. This
  genotypic data has actually been read in from a call to
  \code{wgaim.read.cross} and is therefore a list with components of the
  a returned object from \code{wgaim.read.cross}. See the aformentioned
  function for more details. 
  }
\source{
 
}
\examples{

data(raccas, package = "wgaim")
link.map(raccas, cex = 0.5)

}
\keyword{datasets}
