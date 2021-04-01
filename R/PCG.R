PCG <- function(data, clusters=NULL, alpha=.05, fineTuning=F){

  ### Step 1: Get clusters====
  cat("Step 1: Get clusters\n")
  ## Polychoric correlation
  #require(psych)
  #require(qgraph)
  Corr <- tryCatch(cor_auto(data, forcePD=T), error=function(e) {
          tryCatch(mixedCor(data)$rho, error=function(e) {
            cor(data, method="spearman")
          })
  })
  n <- nrow(data)
  if(length(clusters) == ncol(data)) {
      node_set <- clusters
    } else if (clusters=="ega") {
      #require(EGAnet)
      Dims     <- EGA(data=Corr, n=n, plot.EGA=F)
      node_set <- Dims$wc
  } else if (clusters=="cgmm") {
      #require(mclust)
      HCA      <- Mclust(Corr)
      node_set <- HCA$classification
  } else {"No clusters or estimation method identified!"}

  ### Step 2: Group correlations====
  cat("Step 2: Group correlations\n")
  ## Average correlations
  PowerEdges <- groupVS(Corr,node_set)
  PowerEdges <- PowerEdges[!is.na(rowMeans(PowerEdges, na.rm=T)),
                           !is.na(colMeans(PowerEdges, na.rm=T))]
  V <- rownames(PowerEdges) <- colnames(PowerEdges)
  V <- as.character(V)
  Ndata <- data.frame(apply(data, 2, as.numeric))
  colnames(Ndata) <- colnames(data)

  ### Step 3: Structure Learning====
  cat("Step 3: Structure Learning\n")
  ## IS1: Constraint-based structure learning algorithm
  #require(pcalg)
  skel <- skeleton(suffStat=list(C=PowerEdges, n=n),
                   indepTest=gaussCItest, alpha=alpha,
                   method="stable.fast", labels=V)
  SK <- as(skel,'amat') == 1
  fit <- pc(suffStat=list(C=PowerEdges, n=n), fixedGaps=(SK==F),
            indepTest=gaussCItest, alpha=alpha,
            labels=V, skel.method="stable.fast")
  PC <- as(fit,'amat') == 1

  ## IS2: Get full graph
  # Directed graph
  SKT <- Corr; PCT <- Corr
  for (i in 1:ncol(Corr)) {
    for (j in 1:ncol(Corr)) {
      SKT[i,j] <- tryCatch(SK[node_set[i],node_set[j]], error=function(e) NA)
    }
  }
  for (i in 1:ncol(Corr)) {
    for (j in 1:ncol(Corr)) {
      PCT[i,j] <- tryCatch(PC[node_set[i],node_set[j]], error=function(e) NA)
    }
  }
  # Undirected graph
  UG <- matrix(NA, nrow=length(node_set),ncol=length(node_set))
  colnames(UG) <- row.names(UG) <- colnames(data)
  for(r in 1:length(node_set)) {
    for(c in 1:length(node_set))
      UG[r,c] <- if(node_set[r] == node_set[c]) {1} else {0}
  }

  if (fineTuning == T) {
    ### Step 4: Three types of fine tuning for implied Chain Graphs
    cat("Step 4: Fine-tuning\n")
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

    #require(bnlearn)
    fitp <- pc.stable(Ndata, blacklist = Blist)
    fith <- hc(Ndata, blacklist = Blist)
    fitm <- mmhc(Ndata, blacklist = Blist)
    PCfu <- t(amat(fitp))
    HCfu <- t(amat(fith))
    MMHC <- t(amat(fitm))

    Result <- list("PCG" = t(PC),
                   "clusters"= node_set,
                   "powerEdges" = PowerEdges,
                   "CG_U"=t(UG),
                   "CG_D"=t(PCT),
                   "CG_PC"=PCfu,
                   "CG_HC"=HCfu,
                   "CG_MMHC"=MMHC)
  } else if (fineTuning == F) {
    Result <- list("PCG" = t(PC),
                   "clusters"= node_set,
                   "powerEdges" = PowerEdges,
                   "CG_U"=t(UG),
                   "CG_D"=t(PCT))
  } else stop("Unknow request for fine-tuning.")

  return(Result)
}
