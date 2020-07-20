
# Rt_fun <- function(theta,x){
#   # theta[1]/(1+exp(theta[3]*(x-theta[2])))
#   exp(log(theta[1])-b*x)
# }
# plot(seq(0,1,.01),Rt_fun(theta = c(R0,b),x = seq(0,1,.01)))

# vectorised version
Rt_fun <- function(theta,x){
  R0 <- rep(1,nrow(x)) %*% t(theta[(1-1)*N_geo+ (1:N_geo)])
  # K <- rep(1,nrow(x)) %*% t(theta[(3-1)*N_geo+ (1:N_geo)])
  # X0 <- rep(1,nrow(x)) %*% t(theta[(2-1)*N_geo+ (1:N_geo)])
  B <- rep(1,nrow(x)) %*% t(theta[(2-1)*N_geo+ (1:N_geo)])
  
  res <- exp(log(R0)-B*(1-x))  #R0/(1+exp(K*(x-X0)))
  return(res)
}