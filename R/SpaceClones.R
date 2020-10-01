
SpaceClones3 <- function(X, D, fc_cutoff = log(3/2,base=2),
                         xi_threshold = 10, min_clique_size = 10) {
  UpRegMat <- matrix(0, nrow = ncol(X), ncol = nrow(X))

  rownames(UpRegMat) <- colnames(X)

  for(i in 1:ncol(X)) {
    data <- X[,i]
    it <- which(data > fc_cutoff)

    z_scores <- sapply(data, function(x) {
      (x-mean(data))/sd(data)
    })

    moran_stat <- sapply(c(1:nrow(X)), function(x) {
      z_scores[x]*sum(D[x,]*z_scores)
    })

    myp <- c()

    for(j in it) {
      pval <- perms(D[j,], z_scores, moran_stat[j], j, 200)
      myp <- c(myp, pval)
    }

    #myp <- p.adjust(myp, method = "fdr")
    myp[myp < 0.05] <- 1
    myp[myp < 1] <- 0
    UpRegMat[i,it] <- myp
  }


  cnts <- apply(UpRegMat, 1, function(x) length(which(x != 0)))
  UpRegMat <- UpRegMat[cnts >= xi_threshold, ]

  my.data <- as.matrix(UpRegMat)
  rs <- rowSums(my.data)
  my.data <- my.data[order(rs, decreasing = TRUE),]
  UpRegMat <- UpRegMat[order(rs, decreasing = TRUE),]

  Edgelist <- matrix(nrow = 0, ncol = 2)

  my.binmat <- as.data.frame(UpRegMat)
  rs <- rowSums(my.binmat)
  my.binmat <- my.binmat[order(rs,decreasing = TRUE),]

  pvals <- c()

  if(nrow(my.binmat) == 0) {
    stop("No autocorrelated genes found. No subclones detected.")
  }

  for(i in 1:nrow(my.binmat)) {
    mp <- sum(my.binmat[i,])
    S1 <- which(my.binmat[i,] == 1)
    for(j in i:nrow(my.binmat)) {
      if(i == j) {next}
      #Test if j is a subset of i

      R <- ncol(my.binmat)
      r <- sum(my.binmat[j,])

      k <- length(intersect(which(my.binmat[i,] == 1), which(my.binmat[j,] == 1)))

      pval <- phyper(k, mp, R-mp, r, lower.tail = FALSE)
      pvals <- c(pvals, pval)
      #if(pval > 1 - 10^{-5}) {
      #Edgelist <- rbind(Edgelist, c(i, j))
      #}

      #S2 <- which(my.binmat[j,] == 1)
      #if(sum(S2 %in% S1) >= 0.75*(length(S2))) {
      #  Edgelist <- rbind(Edgelist, c(i,j))
      #}
    }
  }

  if(length(pvals) > 0) {
    pvals <- p.adjust(pvals, method = "fdr")
  }

  cntr <- 1
  for(i in 1:nrow(my.binmat)) {
    mp <- sum(my.binmat[i,])
    S1 <- which(my.binmat[i,] == 1)
    for(j in i:nrow(my.binmat)) {
      if(i == j) {next}
      #Test if j is a subset of i

      if(pvals[cntr] < 0.05) {
        Edgelist <- rbind(Edgelist, c(i,j))
      }

      cntr <- cntr + 1
    }
  }

  G <- graph_from_edgelist(Edgelist, directed = FALSE)
  Cliques <- max_cliques(G)
  cl <- sapply(Cliques, function(x) length(x))


  subclones <- list()
  for(i in 1:nrow(my.data)) {
    l <- sapply(Cliques, function(x) {
      if(i %in% x) {
        return(x)
      }
      else{
        return(NULL)
      }
    })

    lengths <- sapply(l, function(x) length(x))
    ix <- which.max(lengths)
    subclones[[i]] <- l[[ix]]
  }

  lengths <- sapply(subclones, function(x) length(x))

  M <- matrix(0, nrow = nrow(my.data), ncol = nrow(my.data))

  if(nrow(M) == 0) {
    stop("No clique above min_clique_size. No subclones detected.")
  }

  for(i in 1:nrow(M)) {
    M[i,subclones[[i]] ] <- 1
  }

  rownames(M) <- rownames(UpRegMat)
  colnames(M) <- rownames(UpRegMat)
  M <- M[rowSums(M) >= min_clique_size,rowSums(M) >= min_clique_size]
  my.binmat2 <- BinaryMatrix(t(M))
  my.binmat2 <- threshLGF(my.binmat2, cutoff=0.3)

  out <- Mercator(my.binmat2, metric = "binary", method = "tsne", K= my.binmat2@reaper@nGroups, perplexity = 10)

  plot(out)

  return(out)
}

perms <- function(spatial_weights, z_scores, local_moran, i, iters) {
  res <- sapply(1:iters, function(x) {
    newz_scores <- z_scores
    newz_scores[-i] <- sample(z_scores[-i], length(z_scores)-1, replace = FALSE)
    return(z_scores[i]*sum(spatial_weights*newz_scores))
  })

  return(length(which(local_moran <= res))/iters)
}
