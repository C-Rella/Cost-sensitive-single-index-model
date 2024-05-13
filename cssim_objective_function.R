
### --- Optimization parameters --- ###
# As introduced in Klein & Spady (1993): # 0 < b < c; a > 2b + 2c > 0; For e>0 -> a=e, b=e/5, c=e/4
e_si <- 20
e_si_tilde <- e_si*.8
a_si <- e_si
b_si <- a_si/5
c_si <- a_si/4 

### --- Objective function --- ###
cssim_of <- function(theta, bw = length(y)**(-1/5), y, x, w, lambda, a_cost = .1, b_cost = 10, c_cost = 1) {
  k <- ncol(x) 
  n <- nrow(x)
  
  # --- Score for a given theta --- #
  lp <- as.matrix(x)%*%theta
  
  # --- Link function given a score --- #
  g_1 <- npksum(txdat = lp[y == 1], exdat = lp, bws = bw)$ksum 
  g_0 <- npksum(txdat = lp[y == 0], exdat = lp, bws = bw)$ksum 
  p <-  as.vector(g_1/(g_0+g_1))
  
  ### --- Klein & Spady complete proposal --- ###
  # K&S recommend avoid d_i terms to avoid computational errors
  # tau_0 <- ( 1+exp( (bwsi^(e_si_tilde/5)-g_0)/bwsi^(e_si_tilde/4)) )^(-1)
  # tau_1 <- ( 1+exp( (bwsi^(e_si_tilde/5)-g_1)/bwsi^(e_si_tilde/4)) )^(-1)
  # tau <- tau_0*tau_1 
  # z_0 <- (bwsi^b_si - g_0)/bwsi^c_si
  # z_1 <- (bwsi^b_si - g_1)/bwsi^c_si
  # d_0 <- bwsi^a_si*(exp(z_0)/(1+exp(z_0)))
  # d_1 <- bwsi^a_si*(exp(z_1)/(1+exp(z_1)))
  # p <-  (g_1 + d_1)/(g_0+g_1+d_1+d_0)
  
  xi <- -c_cost*y*w + (1-y)*w*a_cost + b_cost # loss function weights
  
  # --- Average expected cost (AEC) --- #
  loss <- sum(c_cost*y*w + p*xi) / n + lambda*sum(abs(theta))  # CS proposal
  return( -loss ) 
  # sum( y*log(p^2) + (1-y)*log((1-p)^2) ) / n # Klein & SPady (1993 proposal)
}

