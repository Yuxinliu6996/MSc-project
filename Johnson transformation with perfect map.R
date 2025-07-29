
Perfect_map_johnsonsu <- function(log_pi,xi,lambda,delta,gamma,y_curr = rep(0,length(xi)),h = 0.9,target_a = 0.23){
  x_from_y <- function(y) sinh((y - gamma)/delta)*lambda + xi 
  
  # Define the transformed density p_Y(y)
  log_piy <- function(y) {
    x <- x_from_y(y)
    logf <- log_pi(x) + sum(log(1 + ((x - xi) / lambda)^2 ) /2)
    if(!is.finite(logf)) logf <- -1e30
    return(logf)}
  
  s <- Adaptive_RWM(log_piy,nits,h,y_curr)
  samples_y = as.matrix(s$x_store)
  samples_x = t(apply(samples_y,1,x_from_y))
  if(dim(samples_x)[1]<dim(samples_x)[2]) samples_x <- t(samples_x)
  return(list(samples_y = samples_y,samples_x = samples_x,step_size = s$step_size))
  
}
