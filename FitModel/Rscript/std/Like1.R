#' get likelihood  
#'
#' get Poisson likelihood  
#' internal to the MCMC
#' 
#' @param lambda: 'force of infection' matrix (incidence weighted by serial interval),
#'                  column number of days, row: number of locations   
#'                  
#' @param I matrix of observed incidence, same dimension as lambda
#' 
#' @param R0 vector of reproduction numbers per locations
#'
#' @details  L log likelihood
#' @export
#' 

Like1 <- function(theta){
  #mobility
  R_daily <- Rt_fun(theta = theta, x = M )
  Rt <- H %*% R_daily

  # logL <- colSums(dpois(x = mD, lambda = Rt*Ot, log = TRUE) ,na.rm=TRUE)
  # logL <- colSums(dnbinom(x = mD, mu = Rt*Ot, size =  theta[2*N_geo+1], log = TRUE) ,na.rm=TRUE)
  logL <- colSums(dnbinom(x = mD, mu = Rt*Ot, size = sqrt(Rt*Ot)* (theta[2*N_geo+1]), log = TRUE) ,na.rm=TRUE)
  # logL <- (dpois(x = mD, lambda = Rt*Ot, log = TRUE) )
  
  return(logL)
}
