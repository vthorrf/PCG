\name{amat2lavaan}
\alias{amat2lavaan}
\title{Adjacency matrix to lavaan syntax}
\description{
Creates `lavaan` syntax for an adjacency matrix estimated from a causal discovery algorithm.
}
\usage{
amat2lavaan(amat)
}
\arguments{
   \item{amat}{The estimated adjacency matrix from a causal discovery algorithm.}
}
\value{
The specification of a linear model using the `lavaan` syntax.
}
\examples{
### Use Big Five data
require(psych)
data("bfi")
data <- bfi[complete.cases(bfi),1:25]

### Run the analysis
fit <- PCG(data, clusters="ega", alpha=.05, fineTuning=F)

### Extract the results
pcg <- fit$PCG
ug  <- fit$CG_U
dg  <- fit$CG_D

### Generate the lavaan code
amat2lavaan(pcg)
}
