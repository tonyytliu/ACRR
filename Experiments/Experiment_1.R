library(pcalg)
library(graph)

p <- 8
n <- 5000

right <- 0
wrong <- 0
right_0.01 <- 0
wrong_0.01 <- 0
right_0.05 <- 0
wrong_0.05 <- 0
right_0.1 <- 0
wrong_0.1 <- 0

res_vec <- c()

for (k in 1:500){
 # set.seed(k)
  gGtrue <- randomDAG(p, prob = 0.5)
  gmG  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
  gmG8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)
  
  
  
  suffStat <- list(C = cor(gmG8$x), n = nrow(gmG8$x))
  test <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
  
  ITRScor_max <- list(c(0,0,0))
  ITRScor_select <- 0
  correct_select <- 0
  resMax <- NULL
  res_matrix <- NULL
  true_matrix <- as(gmG8$g, "matrix")
  alpha_select <- 0.0
  correct_max <- 0
  correct_rate <- 0
  correct_select <- 0
  correct_max_ITRScor <- 0
  correct_max_alpha <- 0
  correct_0.01 <- 0
  correct_0.05 <- 0
  correct_0.1 <- 0
  correct_min <- p*p
  for (alpha in test){
    re_lst <- pc_mod(suffStat, indepTest = gaussCItest, p = ncol(gmG8$x),
                     alpha = alpha, verbose = FALSE, u2pd = "rand") 
    res <- re_lst$res
    res_matrix <- t(as(res, "matrix"))
    ITRScor <- re_lst$ITRScor
    if (is.nan(ITRScor)) {ITRScor = 0}
    correct <- 0
    for (i in 1:p) {
      for (j in i:p){
        if (true_matrix[i, j] > 0){true_matrix[i, j] = 1}
        if (true_matrix[j, i] > 0){true_matrix[j, i] = 1}
        if ((res_matrix[i, j] == true_matrix[i, j]) && (res_matrix[j, i] == true_matrix[j, i])){
          correct <- correct + 1
        }
      }
    }
    
    if (alpha == 0.01) {
      correct_0.01 <- correct}
    if (alpha == 0.05) {
      correct_0.05 <- correct}
    if (alpha == 0.1) {
      correct_0.1 <- correct}
    
    if (correct > correct_max){
      correct_max <- correct
      correct_max_ITRScor <- ITRScor
      correct_max_alpha <- alpha
    }
    if (correct < correct_min){
      correct_min <- correct
    }
    ITRScor_max <- ITRScor_max[order(sapply(ITRScor_max, function(x) x[1], simplify=TRUE), decreasing=FALSE)]
    if (ITRScor > ITRScor_max[[1]][1]){
      ITRScor_max[[1]][1] <- ITRScor
      resMax <- res
      ITRScor_max[[1]][3] <- alpha
      ITRScor_max[[1]][2] <- correct}
  }
  
  ITRScor_max <- ITRScor_max[order(sapply(ITRScor_max, function(x) x[1], simplify=TRUE), decreasing=TRUE)]
  
  res_vec <- c(res_vec, (correct_max - ITRScor_max[[1]][2]))
  res_vec <- c(res_vec, (correct_max - correct_0.01))
  res_vec <- c(res_vec, (correct_max - correct_0.05))
  res_vec <- c(res_vec, (correct_max - correct_0.1))
  
  
  print(k)
}

dim(res_vec) <- c(4, k)
res_vec <- t(res_vec)
write.csv(res_vec, "res_vec.csv")
