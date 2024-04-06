\name{friedman.delgiudice.test}
\alias{friedman.delgiudice.test}
\title{Friedman Del Giudice Test}
\description{
Performs the Friedman Del Giudice test with replications for comparing multiple, replicated treatments across blocks.
}
\usage{
friedman.delgiudice.test(y, treatments, blocks, replications = 1, alpha = 0.05)
}
\arguments{
  \item{y}{A numeric vector of data values in order of block first and treatment second.}
  \item{treatments}{The number of treatments.}
  \item{blocks}{The number of blocks.}
  \item{replications}{The number of replications (default is 1).}
  \item{alpha}{The significance level for the test (default is 0.05).}
  \item{...}{Additional arguments to be passed to the function.}
}
\details{
The Friedman Del Giudice test is a non-parametric test that extends the Friedman test by allowing for replication within each block-treatment combination.
}
\value{
A list with the following components:
\item{statistic}{The test statistic.}
\item{parameter}{The degrees of freedom for the test.}
\item{p.value}{The p-value for the test.}
\item{method}{A character string indicating the name of the test.}
\item{data.name}{A character string giving the name of the data argument.}
\item{conclusion}{A character string indicating whether to "reject" or "fail to reject" the null hypothesis based on the p-value and the specified significance level.}
}
\references{
Spotify & Sleep
}
\author{
Luke Del Giudice
}
\seealso{
\code{\link{friedman.test}} for the Friedman test without replications.
}
\examples{
# Generate example data
y <- c(3, 4, 5, 2, 1, 6, 4, 3, 5, 1, 2, 6, 5, 4, 3, 2, 1, 6)
treatments <- 3
blocks <- 2
replications <- 2

# Perform the Friedman Del Giudice test
result <- friedman.delgiudice.test(y, treatments, blocks, replications)
print(result)
}
\keyword{htest}
