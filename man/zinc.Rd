\name{zinc}
\alias{zinc}
\docType{data}
\title{Zinc data}
\description{Zinc concentration data of a Double Haploid wheat population.
}
\usage{data(zinc)}
\format{A data frame with 200 observations on the following 5 variables
  \describe{
    \item{id}{The identification for the Double Haploid or Control line}
    \item{Block}{A blocking variable for the experiment}
    \item{Type}{The type of wheat variety (Double Haploid, Cascades,
      Rac875-2)}
    \item{shoot}{The shoot length of each plant}
    \item{znconc}{Zinc concentration level of each plant} 
    }
 }
\source{ 
}
\examples{

data(zinc, package = "wgaim")
plot(zinc$shoot, zinc$znconc)

}
\keyword{datasets}
