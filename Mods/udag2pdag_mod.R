#Modified from the package pcalg
#Reference https://CRAN.R-project.org/package=pcalg 

udag2pdag_mod <- function (gInput, verbose = FALSE, rScor, pScor, rScor_count) 
{
  res <- gInput
  if (numEdges(gInput@graph) > 0) {
    g <- as(gInput@graph, "matrix")
    #print(g)
    p <- as.numeric(dim(g)[1])
    pdag <- g
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x)
      for (z in allZ) {
        if (g[x, z] == 0 && !(y %in% gInput@sepset[[x]][[z]] || 
                              y %in% gInput@sepset[[z]][[x]])) {
          if (verbose) {
            cat("\n", x, "->", y, "<-", z, "\n")
            cat("Sxz=", gInput@sepset[[z]][[x]], "Szx=", 
                gInput@sepset[[x]][[z]])
            cat("pScor=", pScor[x, y])
          }
          pdag[x, y] <- pdag[z, y] <- 1
          pdag[y, x] <- pdag[y, z] <- 0
          rScor <- rScor + (pScor[x, y] + pScor[y, z])/2 #add pScor used
          pScor[x, y] <- 0
          pScor[y, z] <- 0
          rScor_count <- rScor_count + 1 #add count
        }
      }
    }
    #print(pdag)
    res2 <- pdag2dag(as(pdag, "graphNEL"))
    if (res2$success) {
      old_pdag <- matrix(0, p, p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE)
        #print(ind)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[b, ] == 1 & pdag[, b] == 
                           1) & (pdag[a, ] == 0 & pdag[, a] == 0))
          if (length(indC) > 0) {
            pdag[b, indC] <- 1
            pdag[indC, b] <- 0
            rScor <- rScor + pScor[b, indC] #add used pScor to rScor
            pScor[b, indC] <- 0
            rScor_count <- rScor_count + 1 #add count
            if (verbose) 
              cat("\nRule 1:", a, "->", b, " and ", b, 
                  "-", indC, " where ", a, " and ", indC, 
                  " not connected: ", b, "->", indC, "\n")
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           0) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) > 0) {
            pdag[a, b] <- 1
            pdag[b, a] <- 0
            rScor <- rScor + pScor[a,b]#add used pScor to rScor
            pScor[a,b] <- 0
            rScor_count <- rScor_count + 1 #add count
            if (verbose) 
              cat("\nRule 2: Kette ", a, "->", indC, 
                  "->", b, ":", a, "->", b, "\n")
          }
        }
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pdag[a, ] == 1 & pdag[, a] == 
                           1) & (pdag[, b] == 1 & pdag[b, ] == 0))
          if (length(indC) >= 2) {
            g2 <- pdag[indC, indC]
            if (length(g2) <= 1) {
              g2 <- 0
            }
            else {
              diag(g2) <- rep(1, length(indC))
            }
            if (any(g2 == 0)) {
              pdag[a, b] <- 1
              pdag[b, a] <- 0
              rScor <- rScor + pScor[a,b]#add used pScor to rScor
              pScor[a,b] <- 0
              rScor_count <- rScor_count + 1 #add count
              if (verbose) 
                cat("\nRule 3:", a, "->", b, "\n")
            }
          }
        }
      }
      res@graph <- as(pdag, "graphNEL")
    }
    else {
      res@graph <- res2$graph
      res@graph <- dag2cpdag(res@graph)
    }
  }
  #print(as(res@graph, "matrix"))
  re_lst <- list(res = res, ITRScor = rScor/rScor_count)
  return(re_lst)
}