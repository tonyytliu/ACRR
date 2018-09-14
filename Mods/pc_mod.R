#Modified from the package pcalg
#Reference https://CRAN.R-project.org/package=pcalg 

pc_mod <- function (suffStat, indepTest, alpha, labels, p, fixedGaps = NULL, 
                    fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = c("relaxed", 
                                                                              "rand", "retry"), skel.method = c("stable", "original", 
                                                                                                                "stable.fast"), conservative = FALSE, maj.rule = FALSE, 
                    solve.confl = FALSE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule) 
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl) 
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule) 
    stop("Choose either conservative PC or majority rule PC!")
  skel_mod<- skeleton_mod(suffStat, indepTest, alpha, labels = labels, 
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges, 
                   NAdelete = NAdelete, m.max = m.max, numCores = numCores, 
                   verbose = verbose)
  skel <- skel_mod$pc_alg
  #record scores from skel_mod
  rScor <- skel_mod$rScor
  pScor <- skel_mod$pScor
  rScor_count <- skel_mod$rScor_count
  #print(as(skel@graph, "matrix"))
  skel@call <- cl
  if (!conservative && !maj.rule) {
    #change rand to record score
    switch(u2pd, rand = udag2pdag_mod(skel, verbose = verbose, rScor = rScor, pScor = pScor, rScor_count = rScor_count), 
           retry = udag2pdagSpecial(skel, verbose = verbose)$pcObj, 
           relaxed = udag2pdagRelaxed(skel, verbose = verbose, 
                                      solve.confl = solve.confl))
  }
  else {
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha, 
                          version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
    udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl, 
                     solve.confl = solve.confl)
  }
}