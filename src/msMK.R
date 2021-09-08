# student.R

## usethis namespace: start
#' Imports
#' @importFrom Rcpp sourceCpp
#' @useDynLib msMK, .registration = TRUE
#' @export bintree
## usethis namespace: end
"_PACKAGE"

Rcpp::loadModule(module = "msMK", what = TRUE)
