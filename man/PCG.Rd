\name{PCG}
\alias{PCG}
\title{Power Chain Graph Structure Learning}
\description{
Learns a power chain graph following Franco et al. (2022).
}
\usage{
PCG(data, cor="auto", clusters=NULL, alpha=.05, fineTuning=F)
}
\arguments{
   \item{data}{A data frame of the variables to be included in the PCG.}
   \item{cor}{The correlation method to be used. Defaults to 'auto'. Can also be "pearson", "spearman", "kendall", or "poly".}
   \item{clusters}{A integer vector of expected classes for each variable in data. If an exploratory approach is preferred, use "ega" for learning the cluster by the EGA method or "cgmm" learning the cluster by the CGMM method.}
   \item{alpha}{The significance level for testing the directional relations with the PC-stable algorithm.}
   \item{fineTuning}{Boolean indicating if fine-tuning should be performed. Defaults to FALSE.}
}
\value{
A list containing the following component:
  \item{PCG}{The estimated PCG.}
  \item{clusters}{The integer vector of expected classes for each variable in data.}
  \item{CG_U}{The undirected part of the chain graph induced by the PCG.}
  \item{CG_D}{The directed part of the chain graph induced by the PCG.}
  \item{CG_PC}{The directed part of the chain graph implied by the PCG, fine-tuned with the PC-stable algorithm.}
  \item{CG_HC}{The directed part of the chain graph implied by the PCG, fine-tuned with the hill-climbing algorithm.}
  \item{CG_MMCC}{The directed part of the chain graph implied by the PCG, fine-tuned with the Min-Max hill-climbing algorithm.}
}
\references{
Franco, V. R., Barros, G., Wiberg, M., & Laros, J. A. (2022). Chain graph reduction into power chain graphs. Quantitative and Computational Methods in Behavioral Sciences, 2(1). https://doi.org/10.5964/qcmb.8383
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

### Visualize
require(qgraph)
require(EGAnet)
par(mfrow=c(2,1))
qgraph::qgraph(pcg)
qgraph::qgraphMixed(ug * qgraph::cor_auto(data, forcePD=T),
                    dg * EGAnet::EGA(data, plot=F)$network,
                    groups = as.factor(fit$clusters),
                    labels = colnames(data))
}
