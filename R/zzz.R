.First.lib <- function(lib, pkg){
  if(!exists("asreml.Rsys", envir = .GlobalEnv)) {
    message("ASReml-R needs to be installed before this package can be used.

Please visit http://www.vsni.co.uk/products/asreml/ for more information.\n")
      }
}
  

