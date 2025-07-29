
# Define the fitting function in one dimension using MLE
fit_oned_johnsonsu_ML <- function(samples){
  # Define a negative log-likelihood function of jognson su distribution 
  loglik_johnson_su <- function(par, x) {
    gamma <- par[1]
    delta <- par[2]
    xi    <- par[3]
    lambda<- par[4]
    
    pi <- 3.1415926535
    if (delta <= 0 || lambda <= 0) return(Inf)
    
    y <- (x - xi) / lambda
    z <- gamma + delta * asinh(y)
    
    log_pdf <- log(delta) - log(lambda) - 0.5 * log(2*pi) - 0.5 * z^2 + log(1 / sqrt(1 + y^2))
    log_pdf <- pmax(log_pdf, -1e30)
    return(-sum(log_pdf))  # pay attention it is negative 
  }
  start <- c(0, 1, 0, 1)
  
  # calculate the maximum likelihood estimator
  result_ml <- optim(par=start, fn=loglik_johnson_su, x=samples, method="L-BFGS-B", 
                     lower=c(-Inf, 1e-6, -Inf, 1e-6))
  # extrat it from the result
  gamma <- result_ml$par [1]
  delta <- result_ml$par [2]
  xi    <- result_ml$par [3]
  lambda<- result_ml$par [4]
  
  # return the parameters
  c(xi, lambda, gamma, delta)
}
# Find the approximation in high dimension using component-wise method
fit_highdim_johnsonsu_ML <- function(X) {
  if(!is.matrix(X)){
    return("the input must be a matrix")
  }
  # Find the Johnson approximation in each dimension
  r <- as.data.frame(t(apply(X,2,fit_oned_johnsonsu_ML)))
  colnames(r) <- c("xi","lambda","gamma","delta")
  return(r)
}
transformation_highdim_map_johnsonsu <- function(xi,lambda,gamma,delta,nits,log_pi,x_curr = rep(0,length(xi)), target_a= 0.44){
  d <- length(xi)
  # Define inverse map: x = F⁻¹(Φ(y))
  
  x_from_y <- function(y) sinh((y - gamma)/delta)*lambda + xi 
  
  # Define the transformed density p_Y(y)
  log_p_Y <- function(y) {
    x <- x_from_y(y)
    logf <- log_pi(x) + sum(log(1 + ((x - xi) / lambda)^2 ) /2)
    logf <- max(logf,-1e30)
    if(is.nan(logf)) logf <- -1e30
    return(logf)
  }
  
  # Sampling from log_p_Y
  chain<- Adaptive_RWM(log_p_Y,nits = nits,h = 0.5,x_curr = x_curr,target_a = target_a)
  
  # Get the samples of y
  samples_y <- chain$x_store
  
  # transform y back
  samples_x<- t(apply(samples_y,1,x_from_y)) 
  if(dim(samples_x)[1]<dim(samples_x)[2]) samples_x = t(samples_x)
  
  return(list(samples_x = samples_x,samples_y = samples_y,
              a_rate = chain$a_rate,step_size = chain$step_size))
}
Sampling_trsanformation_johnsonsu_MLE <- function(samples,log_pi,nits,method = "MLE",x_curr = rep(0,ncol(samples)), target_a = 0.23){
  if(method == "MLE") johnsonsu_approx <- fit_highdim_johnsonsu_ML(samples)
  if(method == "Moment") johnsonsu_approx <- fit_highdim_johnsonsu_moment(samples)
  xi <- johnsonsu_approx$xi
  lambda <- johnsonsu_approx$lambda
  gamma <- johnsonsu_approx$gamma
  delta <- johnsonsu_approx$delta
  result <- transformation_highdim_map_johnsonsu(xi,lambda,gamma,delta,nits,log_pi,x_curr,target_a)
  return(list(xi = xi,lambda = lambda,gamma = gamma,delta = delta,
              samples_x = result$samples_x,
              samples_y = result$samples_y,
              a_rate = result$a_rate,step_size = result$step_size))
}
