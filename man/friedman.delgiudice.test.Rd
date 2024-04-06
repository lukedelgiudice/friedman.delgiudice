\name{friedman.delgiudice.test}
\alias{friedman.delgiudice.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
place holder
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
friedman.delgiudice.test(y, treatments, blocks, replications = 0, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{treatments}{
%%     ~~Describe \code{treatments} here~~
}
  \item{blocks}{
%%     ~~Describe \code{blocks} here~~
}
  \item{replications}{
%%     ~~Describe \code{replications} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, treatments, blocks, replications = 0, alpha = 0.05,
    ...)
{
    t <- treatments
    b <- blocks
    r <- replications + 1
    if (anyNA(y)) {
        stop("NA's are not allowed in 'y'")
    }
    if (length(y) <= 0) {
        stop("'y' must have data")
    }
    if (any(c(r, t, b)) <= 0L) {
        stop("'treatments' 'blocks', and 'replications' must be valid inputs")
    }
    if (length(y) != (t * b * r)) {
        stop("dimensions of the given data and arguments don't match")
    }
    blockdf <- data.frame(y, ncol = (t * r), byrow = T)
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
    test_statistic <- (12 * b)/((t * r) * ((t * r) + 1)) * sum((average_treatment_ranks -
        ((t + 1)/2))^2)
    p_value <- pchisq(test_statistic, t - 1, lower.tail = FALSE)
    result <- list(statistic = test_statistic, parameter = t -
        1, p.value = p_value, method = "Friedman Del Giudice Test",
        data.name = deparse1(substitute(y)), conclusion = ifelse(p_value <
            alpha, "reject", "fail to reject"))
    class(result) <- "htest"
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.