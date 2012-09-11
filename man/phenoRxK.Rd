\name{phenoRxK}
\alias{phenoRxK}
\docType{data}
\title{Phenotypic RAC875 x Kukri trial data frame}
\description{Phenotype data arising from a field trial of a
  Double Haploid population involving a crossing of the wheat varities
  RAC875 and Kukri}
\usage{data(phenoRxK)}
\format{This data relates to a field trial conducted in 2007 at the
  Roseworthy Cmapus of the University of Adelaide. The trial consisted
  of 2 replicates of 367 Double Haploid lines from a cross between wheat
  varieties RAC875 and Kukri. The DH lines, the parents
  (RAC875, Kukri) and commercial varieties (Co5693*E002, DRYSDALE, EXCALIBUR,
  KRICHAUFF, WYALKATCHEM, YITPI) were randomly allocated to 756
  plots using a randomized complete bloack design. The trial was laid out
  in a 63 by 12 rectangular array. The data frame consists of 756 rows
  with 8 columns described by
  \describe{
    \item{Genotype:}{A 375 level factor containing a unique
      identification for the wheat varieties involved in the experiment.}
    \item{Type:}{A factor of four levels indexing the wheat varieties
      (Double Haploid, RAC875, Kukri, Other) where Other consist of
      (Co5693*E002, DRYSDALE, EXCALIBUR, KRICHAUFF, WYALKATCHEM, YITPI)}
    \item{id:}{A factor of 183 levels uniquely identifying the wheat
      varieties involved in the experiment.}
    \item{Range:}{A factor of 12 numeric levels indexing the field Range.}
    \item{Row:}{A factor of 63 numeric levels indexing the field Rows.}
    \item{Block:}{A factor of 2 levels indexing the Blocks of the
      experiment}
    \item{yield:}{A numeric vector of yield observations in kg/ha}
    \item{tgw:}{A numeric vector of thousand grain weight observations}
     }
  }
\examples{
data(phenoSxT, package = "wgaim")

}
\keyword{datasets}
