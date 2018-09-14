library(pcalg)
library(graph)

output_vec <- c()

suffStat <- list(C = cor(dataset), n = nrow(dataset), n.min = nrow(dataset),
                 dm = dataset, adaptDF = FALSE) #replace dataset with the dataset you have

test <- c(0.05,0.1) #input the alphas you want to test.

ITRScor_max <- list(c(0,0,0))
resMax <- NULL
res_matrix <- NULL

for (alpha in test){
  #change indepTest to the test correspond to the data set
  re_lst <- pc_mod(suffStat, indepTest = gaussCItest, labels=colnames(dataset),
                   alpha = alpha, verbose = FALSE, u2pd = "rand")
  
  res <- re_lst$res
  res_matrix <- t(as(res, "matrix"))
  ITRScor <- re_lst$ITRScor
  if (is.nan(ITRScor)) {ITRScor = 0}
  
  ITRScor_max <- ITRScor_max[order(sapply(ITRScor_max, function(x) x[1], simplify=TRUE), decreasing=FALSE)]
  if (ITRScor > ITRScor_max[[1]][1]){
    ITRScor_max[[1]][1] <- ITRScor
    resMax <- res
    ITRScor_max[[1]][3] <- alpha}
  cat("\n", alpha, "; ", ITRScor)
  
  output_vec <- c(output_vec, alpha, ITRScor)
}

if (require(Rgraphviz)) {
  plot(resMax, main = "")
}