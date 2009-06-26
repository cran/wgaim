\name{wgaim.asreml}
\alias{wgaim.asreml}
\alias{wgaim}
\title{wgaim method for class "asreml"}
\description{Fits an iterative Whole Genome Average Interval Mapping
  (wgaim) model for QTL detection 
}
\usage{
\method{wgaim}{asreml}(baseModel, parentData, TypeI = 0.05, attempts = 5,
                       trace = TRUE, \ldots)
}
\arguments{
  \item{baseModel}{A model of class "\code{asreml}" usually representing
    a base model with which to build the qtl model.  
  }
  \item{parentData}{A data structure that inherits the class
    "\code{interval}" containing the genotypic data as well as the
    phenotypic data used in the \code{baseModel}. This should contain a
    structure called "\code{full.data}" from appropriately merging
    phenotypic and genotypic data (see \code{read.interval} and \code{wmerge}).     
  }
  \item{TypeI}{The level of significance for detecting a QTL. The default is
  0.05.}
  \item{attempts}{The number of attempts at convergence for the fixed or
    random qtl model. The default is 5.}
  \item{trace}{An automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace machanisms will display a
  message if a QTL is detected.}
  \item{\ldots}{Any other extra arguments to be passed to each of the
    \code{asreml} calls. These may also include \code{asreml.control} arguments.}
}
\details{
The \code{parentData} should contain a data structure "\code{full.data}"
which consists of a genotypic and phenotypic information appropriately merged
(see \code{wmerge}.)   
  
As the Whole Genome Average Interval Mapping approach detects QTL
sequentially it may require many calls to \code{asreml}. For this reason
the function may seem slow and users should be patient.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.

 To avoid lexcial scoping problems the function also places a version of
the "\code{full.data}"  from \code{parentData} into a local data object called
"\code{asdata}". This ensures that the final model may also be
investigated using the methods and features of the \code{asreml}
package.   
}
\value{An object of class "\code{wgaim}" which also inherits the class
  "\code{asreml}" by default. The object returned is actually an \code{asreml}
  object (see \code{asreml.object}) with the addition of components from
  the QTL detection listed below.
  
  \item{QTL}{A list of components from the significant QTL detected
    including a character vector of the significant chromosomes the
    QTL were detected in \code{sig.chr} and the significant intervals
    on each chromosome in \code{sig.int}. For exploratory purposes the
    final random effects qtl model is also returned in the component
    \code{qtlModel}. There are also a number of diagnostic meausres that
    can be found in \code{diag} that are used in conjunction with \code{tr.wgaim}}
}
\references{
 Verbyla, A. P., Cullis, B. R., Thompson, R (2007) The analysis of QTL
 by simultaneous use of the full linkage map. \emph{Theoretical And
    Applied Genetics}, \bold{116}, 95-111.
}
\author{Julian Taylor, Simon Diffey, Ari Verbyla and Brian Cullis}
\seealso{\code{\link{wmerge}}, \code{\link{print.wgaim}}, \code{\link{summary.wgaim}}}
\examples{

\dontrun{
# read in data 

data(zinc, package = "wgaim")
data(raccas, package = "wgaim")

# subset linkage map and merge genotypic with phenotypic

raccasS <- subset(raccas, chr = c("1A1", "2D1", "4D2", "6A1"))
raccasM <- wmerge(raccasS, zinc, by = "id")

# base model

zn.fm <- asreml(znconc ~ Type, random = ~ Block + id, data = zinc)

# find QTL's

zn.qtl <- wgaim(zn.fm, parentData = raccasM, trace = "trace.txt",
na.method.X = "include") 
}

}
\keyword{regression}
