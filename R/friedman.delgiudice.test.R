friedman.delgiudice.test <- function(y, treatments, blocks, replications = 0, alpha = 0.05, ...) {

  t <- treatments
  b <- blocks
  r <- replications + 1

  if (anyNA(y)) {
    stop("NA's are not allowed in 'y'")
  }

  if (length(y) <= 0) {
    stop("'y' must have data")
  }

  if (any(c(r, t, b) <= 0L)) {
    stop("'treatments' 'blocks', and 'replications' must be valid inputs")
  }

  if (length(y) != (t * b * r)) {
    stop("dimensions of the given data and arguments don't match")
  }

  blockdf <- data.frame(y, ncol = (t*r), byrow = T)
  rankdf <- as.data.frame(t(apply(blockdf, MARGIN = 1, rank)))

  average_treatment_ranks <- numeric(t)

  get_indices <- function(subset_num) {
    start_index <- (subset_num - 1) * r + 1
    end_index <- start_index + r - 1
    return(start_index:end_index)
  }

  for (i in 1:t) {
    average_treatment_ranks[i] <- mean(unlist(rankdf[, get_indices(i)]))
  }

  test_statistic <- (12*b) / ((t*r) * ((t*r) + 1)) * sum((average_treatment_ranks - ((t + 1) / 2))^2)
  p_value <- pchisq(test_statistic, t - 1, lower.tail = FALSE)

  result <- list(
    statistic = test_statistic,
    parameter = t - 1,
    p.value = p_value,
    method = "Friedman Del Giudice Test",
    data.name = deparse1(substitute(y)),
    alpha = alpha,
    conclusion = ifelse(p_value < alpha, "reject", "fail to reject")
  )

  class(result) <- "friedman_test"
  return(result)
}

print.friedman_test <- function(x, ...) {
  cat("Friedman Del Giudice Test\n\n")
  cat("data:  ", x$data.name, "\n")
  cat("Friedman chi-squared = ", round(x$statistic, 4), ", df = ", x$parameter, ", p-value = ", format.pval(x$p.value, digits = 4), "\n\n", sep = "")
  cat("Conclusion: ", ifelse(x$conclusion == "reject", "Reject", "Fail to reject"), " the null hypothesis at alpha = ", x$alpha, "\n", sep = "")
}
