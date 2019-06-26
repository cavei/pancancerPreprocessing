#' Impute values
#'
#' Wrapper to impute.knn. ?inpute.knn to see arguments.
#' @param \dots all arguments are passed to impute.knn
#'
#' @importFrom impute impute.knn
#'
#' @export
#'
imputeKNN <- function(...) {
  sink("/dev/null")
  data <- impute::impute.knn(...)
  sink()
  data
}
