\name{Power Chain Graph}
\alias{Power Chain Graph}
\title{Power Chain Graph Structure Learning}
\description{
Learns a power chain graph following Franco, Wiberg & Laros (under review).
}
\usage{
PCG(data, clusters=NULL, alpha=.05)
}
\arguments{
   \item{data}{A data frame of the variables to be included in the PCG.}
   \item{clusters}{A integer vector of expected classes for each variable in data. If an exploratory approach is preferred, use "ega" for learning the cluster by the EGA method or "cgmm" learning the cluster by the CGMM method.}
   \item{alpha}{The significance level for testing the directional relations.}
\value{
A list containing the following component:
  \item{PCG}{The estimate PCG.}
  \item{CG_F}{The chain graph implied by the PCG.}
  \item{CG_PC}{The chain graph implied by the PCG, fine-tuned with the PC-stable algorithm.}
  \item{CG_HC}{The chain graph implied by the PCG, fine-tuned with the hill-climbing algorithm.}
  \item{CG_MMCC}{The chain graph implied by the PCG, fine-tuned with the Min-Max hill-climbing algorithm.}
}
\examples{
### Use Big Five data
require(psych)
data("bfi")
data <- bfi[,1:25]

### Run the analysis
fit <- PCG(data, clusters="ega", alpha=.01)

### Extract the results
pcg <- fit$PCG

### Visualize
require(qgraph)
qgraph(pcg)
}
