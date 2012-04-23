.onAttach <- function(libname, pkgname){
  if(!("asreml" %in% loadedNamespaces())) {
    packageStartupMessage("ASReml-R needs to be installed before this package can be used.

Please visit http://www.vsni.co.uk/products/asreml/ for more information.\n")
      }
}

