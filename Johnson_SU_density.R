
# define the log density
library(coda)
log_dJohnson <- function(x,lambda,delta,xi,gamma){
  z <- (x - xi) / lambda
  term1 <- log(delta) - log(lambda) - log(2 * pi)/2
  term2 <- -0.5 * log(1 + z^2)
  term3 <- -0.5 * (gamma + delta * asinh(z))^2
  log_pdf <- sum(term1 + term2 + term3)
  #print(term1+term2+term3)
  return(max(log_pdf,-1e30))
}

# define the gradient of logpi
d_logJohnson <- function(x,lambda,delta,xi,gamma){
  z <- (x - xi) / lambda
  result <- - 1 / lambda * (z / (1 + z^2) + 1 / sqrt(1 + z^2) * delta *(gamma + delta * asinh(z)))
  return(result)
}
