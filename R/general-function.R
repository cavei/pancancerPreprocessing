#' Expand collapsed "cols" aìin a DataFrame df
#'
#' @param df a DataFrame
#' @param cols the colums that need to be expanded
#' @param sep separator character
#' @param pairs logical; split as pairs
#'
#' @return DataFrame
#' @rdname general-functions
#' @importFrom checkmate assert_class
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#'
#' @export
#'
set_expand <- function(df, cols, sep=";", pairs=TRUE){
  checkmate::assert_class(df,classes = "DataFrame")

  if (is.numeric(cols)) {
    stop("Colnames needed")
  } else {
    unk <- setdiff(cols, colnames(df))
    if (length(unk))
      stop(paste0("Some columns were not found: ", paste(unk, collapse=", ")))
  }

  expandedCols <- lapply(cols, function(c) {
    strsplit(df[[c]], sep)
  })

  expColsList <- lapply(expandedCols, function(x) {
    do.call(c, x)
  })

  names(expColsList) <- cols

  col1list <- lapply(expandedCols[[1]], length)

  if (length(cols) > 1 & pairs) {
    a <- lapply(expandedCols[-c(1)], function(x) {
      collist <- lapply(x, length)
      if (!identical(collist, col1list))
        stop("Pair mode needs the same length when columns are expanded.")
      collist
    })
  }

  mCols <- setdiff(colnames(df), cols)
  mCOlList <- lapply(mCols, function(m) {
    rep(df[[m]], times=col1list)
  })
  names(mCOlList) <- mCols

  DataFrame(mCOlList, expColsList)[,colnames(df)]
}
