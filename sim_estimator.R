library(maxLik)
library(np)

sim_estimator <- function(theta, y, x, x_test) {
  k <- ncol(x) 
  n <- nrow(x)
  n_train <- n
  n_test <- nrow(x_test)
  theta_aux <- c(1, theta[2:k])
  bwsi <- abs(theta[k+1])
  
  lp <- as.matrix(x) %*% theta_aux
  lp_test <- as.matrix(x_test) %*% theta_aux
  g_1 <- npksum(txdat = lp[y == 1], exdat = c(lp, lp_test), bws = bwsi)$ksum 
  g_0 <- npksum(txdat = lp[y == 0], exdat = c(lp, lp_test), bws = bwsi)$ksum 
  p_aux <-  g_1/(g_0+g_1) 
  p <- p_aux[1:n_train]
  p_test <- p_aux[(n_train + 1):(n_train+n_test)]
  
  return( list(p = p, index = lp, p_test = p_test, index_test = lp_test) )
}
