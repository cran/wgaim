\name{news}
\title{Release information for wgaim}
\section{Changes in version 1.4-x}{
\subsection{README}{
  \itemize{
  \item Since this is the first documented NEWS release on the wgaim package
  it contains new features that have been included over several
  versions. In contrast the documented bug fixes are only recent.
 }
}
  \subsection{NEW FEATURES}{
    \itemize{
      \item The package now includes a PDF reference manual that is accessible by
      navigating to the \code{"doc"} directory of the package. This can
      be found on any operating system using the command

      \code{> system.file("doc", package = "wgaim")}

      The reference manual contains WGAIM theory and two thorough examples that show the
      features of the package. It also contains a "casual walk through"
      the package providing the user with a series of 5 steps to a successful wgaim
      analysis.
      \item The package now includes three fully documented phenotypic
      and genotypic data sets for users to explore. Two of these three
      have been used in the manual and scripts that follow the
      examples in the manual are available under the "doc" directory of the package.
      \item The package now provides very efficient whole genome QTL analysis of high
      dimensional genetic marker data. All genetic marker data is passed
      into \code{wgaim.asreml()} through the \code{"intervalObj"}
  argument. Merging of genotypic and phenotypic data occurs within
  \code{wgaim.asreml()}.
      \item \code{wgaim.asreml()} has several new arguments related
  to selection of QTL. The \code{"gen.type"} argument allows the user to
  choose a whole genome marker analysis or whole genome mid-point
  interval analysis from Verbyla et. al (2007). The \code{"method"} argument gives you the choice of placing
  QTL in the fixed part of the linear mixed model as in Verbyla et.al
  (2007) or the random part of model as in Verbyla et. al
  (2012). Finally, the \code{"selection"} argument allows you to choose whether QTL selection
  is based on whole genome interval outlier statistics or a two stage process of
  using chromosome outlier statistics and then interval outlier
  statistics.
     \item A \code{"breakout"} argument is now also provided which allows
  the user to breakout of the forward selection algorithm at any
  stage. The current model along with any calculated QTL components are
  all available for inspection.
  \item All linkage map plotting functions can be subsetted by
  predefined distances. This includes a list of distances as long as the
  number of linkage groups plotted.
}
}
\subsection{BUG FIXES}{
  \itemize{
    \item \code{cross2int()} now accepts R/qtl objects with cross type
    \code{"bc","dh","rilsef"}.
    \item Fixed an issue with the internal function \code{fix.map()} that
    allowed some co-located sets of markers to appear in the final
    reduced linkage map.
    \item Fixed a long standing scoping issue with different versions of ASReml-R.
    \item Fixed an elusive problem that causes wgaim models to increase the size
    of your .RData upon saving. This is actually an inherent problem with
    using model formula within a function a returning the subsequent
    model. There is now a function at the very tail of
    \code{wgaim.asreml()} that quietly destroys the useless environments
    that these formula contain.
    \item Fixed bug that caused \code{wgaim.asreml()} to crash when no QTL
    were found.
    \item Fixed bug that caused \code{summary.wgaim()} to crash when one
    QTL was found using \code{method = "random"}.
  }
}
}


