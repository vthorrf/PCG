groupVS <- function(Corr,blocks) {
  Bs <- length(unique(blocks))
  Fisher <- fisherz(Corr)
  diag(Fisher) <- rep(NA, length(diag(Fisher)))
  Test <- matrix(NA,ncol=Bs,nrow=Bs)
  row.names(Test) <- colnames(Test) <- sapply(1:length(unique(blocks)),
                                              function(g) paste("P",g,sep=""))
  for (i in 1:Bs) {
    for (j in 1:Bs) {
      Test[i,j] <- fisherz2r(mean(Fisher[which(blocks == i),
                                         which(blocks == j)],na.rm=T))
    }
  }
  diag(Test) <- 1
  return(Test)
}
