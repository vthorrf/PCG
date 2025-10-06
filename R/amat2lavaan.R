amat2lavaan <- function(amat) {
  # Control condition
  if(any(c(is.null(rownames(amat)), is.null(colnames(amat))))) {
    stop("The rows and the columns of the adjacency matrix should be named")
  }
  var_names <- colnames(amat)

  # Pre-setting
  n <- ncol(amat)
  model_syntax <- c()

  # Process directed edges (grouped by target)
  for (target in seq_len(n)) {
    predictors <- which(amat[, target] == 1 & amat[target, ] == 0)
    if (length(predictors) > 0) {
      reg <- paste0(
        var_names[target],
        " ~ ",
        paste(var_names[predictors], collapse = " + ")
      )
      model_syntax <- c(model_syntax, reg)
    }
  }

  # Process undirected edges (covariances)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (amat[i, j] == 1 && amat[j, i] == 1) {
        cov <- paste0(var_names[i], " ~~ ", var_names[j])
        model_syntax <- c(model_syntax, cov)
      }
    }
  }

  # Combine all into one string
  lavaan_model <- paste(model_syntax, collapse = "\n")

  return(lavaan_model)
}
