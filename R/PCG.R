PCG <- function(data, clusters=NULL, alpha=.05){

  ### Step 1: Get clusters====
  ## Polychoric correlation
  require(psych)
  Corr <- tryCatch(psych::mixed.cor(data)$rho, error=function(e) {
                   cor(data, method="spearman")
  })
  n <- nrow(data)
  if(length(clusters) == ncol(data)) {
      node_set <- clusters
    } else if (clusters=="ega") {
      require(EGAnet)
      Dims <- EGAnet::EGA(data=data, plot.EGA=F)
      node_set <- Dims$wc
  } else if (clusters=="cgmm") {
      require(mclust)
      HCA    <- mclust::Mclust(Corr)
      node_set <- HCA$classification
  } else {"No cluster values or method identified!"}

  ### Step 2: Group correlations====
  ## Average correlations
  groupVS <- function(Corr,blocks) {
    Bs <- length(unique(blocks))
    Fisher <- psych::fisherz(Corr)
    diag(Fisher) <- rep(NA, length(diag(Fisher)))
    Test <- matrix(NA,ncol=Bs,nrow=Bs)
    row.names(Test) <- colnames(Test) <- sapply(1:length(unique(blocks)),
                                                function(g) paste("P",g,sep=""))
    for (i in 1:Bs) {
      for (j in 1:Bs) {
        Test[i,j] <- psych::fisherz2r(mean(Fisher[which(blocks == i),
                                                  which(blocks == j)],na.rm=T))
      }
    }
    return(Test)
  }
  PowerEdges <- groupVS(Corr,node_set)
  V <- colnames(PowerEdges)
  Ndata <- data.frame(sapply(1:ncol(data), function(g) as.numeric(data[,g])))
  colnames(Ndata) <- colnames(data)

  ### Step 3: Structure Learning====
  ## IS1: Constraint-based structure learning algorithm
  require(pcalg)
  skel <- pcalg::skeleton(suffStat=list(C=PowerEdges, n=n),
                          indepTest=gaussCItest, alpha=alpha,
                          method="stable.fast", labels=V )
  SK <- as(skel,'amat') == 1
  fit <- pcalg::pc(suffStat=list(C=PowerEdges, n=n), fixedGaps=(SK==F),
                   indepTest=gaussCItest, alpha=alpha,
                   labels=V, skel.method="stable.fast" )
  PC <- as(fit,'amat') == 1

  ## IS2: Get full graph
  SKT <- Corr; PCT <- Corr
  for (i in 1:ncol(Corr)) {
    for (j in 1:ncol(Corr)) {
      SKT[i,j] <- SK[node_set[i],node_set[j]]
    }
  }
  for (i in 1:ncol(Corr)) {
    for (j in 1:ncol(Corr)) {
      PCT[i,j] <- PC[node_set[i],node_set[j]]
    }
  }

  ### Step 4: Three types of fine tunning for implied Chain Graphs
  m <- PCT
  Blist1 = data.frame(to=rownames(m)[row(m)[lower.tri(m)]],
                      from=colnames(m)[col(m)[lower.tri(m)]],
                      corr=m[lower.tri(m)])
  Blist1 <- Blist1[,c(2,1,3)]
  Blist2 = data.frame(to=rownames(m)[row(m)[upper.tri(m)]],
                      from=colnames(m)[col(m)[upper.tri(m)]],
                      corr=m[upper.tri(m)])
  Blist2 <- Blist2[,c(2,1,3)]
  Blist <- rbind(Blist1,Blist2)
  Blist <- Blist[Blist$corr == F,-3]

  require(bnlearn)
  fitp <- bnlearn::pc.stable(Ndata, blacklist = Blist)
  fith <- bnlearn::hc(Ndata, blacklist = Blist)
  fitm <- bnlearn::mmhc(Ndata, blacklist = Blist)
  PCfu <- t(bnlearn::amat(fitp))
  HCfu <- t(bnlearn::amat(fith))
  MMHC <- t(bnlearn::amat(fitm))

  Result <- list("PCG" = t(PC),
                 "CG_F"=t(PCT),
                 "CG_PC"=PCfu,
                 "CG_HC"=HCfu,
                 "CG_MMHC"=MMHC)
  return(Result)
}
