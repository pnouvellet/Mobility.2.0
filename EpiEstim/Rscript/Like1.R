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
  #
  Rt <- matrix(rep(theta[1:(N_geo*N_week)],each=7),nrow = N_days, ncol = N_geo, byrow = FALSE)

  # logL <- colSums(dpois(x = mD, lambda = Rt*Ot, log = TRUE) ,na.rm=TRUE)
  logL <- colSums(dnbinom(x = mD, mu = Rt*Ot, size = sqrt(Rt*Ot)* theta[(N_geo*N_week)+1] , log = TRUE) ,na.rm=TRUE)
  # logL <- (dpois(x = mD, lambda = Rt*Ot, log = TRUE) )
  
  return(logL)
}
