#' MCMC iterate
#'
#' run the MCMC to sample posterior of R and initial coniditions at each location 
#' FYI: this is called internally by adapt_tuning
#' 
#' @param incidence the incidence for the time window during which we assume Rt to be constant.  
#'           I is a dataframe, first column are dates then incidence for all locations
#'           nb of row is the size of time widows, dates must be sequential
#' 
#' @param N_geo integer of  numbers of locations
#'                   
#' @param iter integer, the number of iteration for the MCMC
#'
#' @param theta0 vector of inital parameters, here taken from the last MCMC iteration after tuning (save some burn-in)
#'
#' @param s variance of proposal distributions (log-normal) - tuned previously
#' 
#' @param SI Serial interval distribution (see SI_gamma_dist_EpiEstim)
#' 
#' @param mu0: initial conidtions to guaranty that if R=1, then we predict the number of cases in the future will stablise at the mean number of cases observed in the time window
#'              mu0 is also used as the mean of the (exponential) prior for intial conditions estimated
#' 
#' @details  res a list containing 2 matrices: theta: matrix of posterior samples
#'                      and logL: matrix of associated log-likelihood
#' @export
#' 

MCMC_iter <- function(iter,theta0,s){
  
  
  #############################################################################
  # for MCMC
  thetas <- matrix(NA,iter,length(theta0))  # Rt's and initial conditions for each location
  L <- thetas                         # store likelihoods
  
  #############################################################################
  # get the daily 'forces of infections'
  
  
  
  # L1 <- Like1(theta = theta0)
  L1 <- Like1(theta = theta0)
  
  L[1,] <- rep(rep(L1,each = N_week),2)
  thetas[1,] <- theta0       
  
  #############################################################################
  # sampling
  for (i in 2:iter){   
    #print(i)
    for (j in 1:(N_week+1)){
      Ts <- theta0
      if(j == (N_week+1)){
        J <- seq(N_week*N_geo+1,N_week*N_geo*2,by=1)
      }else{
        J <- seq(j,length(theta0)/2,by=N_week)#(j-1)*N_geo+(1:N_geo)
      }
      # propose new parameter j
      if (j == (N_week+1)){
        Ts[J] <- Ts[J]*exp(s[J]*rnorm(1,0,1))
      }else{
        Ts[J] <- Ts[J]*exp(s[J]*rnorm(N_geo,0,1))
      }
      
      f <- which((Ts < prior_theta[,1]) | (Ts > prior_theta[,2]))
      Ts[f] <- theta0[f]
      
      Lint <- Like1(theta = Ts)
      
      
      #get the ratio with previous value of parameter and correct for porposal (and, only for initial conditions, prior distribution)
      # if (j ==1){
      r <- exp(Lint-L1)*Ts[J]/theta0[J]  * exp(-1*(Ts[(N_week*N_geo+1)]-theta0[(N_week*N_geo+1)]))#
      # }else{
      # r <- exp(Lint-L1)
      # }
      # accept or reject
      if (j == (N_week+1)){
        r <- exp(sum(Lint)-sum(L1))*Ts[J[1]]/theta0[J[1]] * exp(-1*(Ts[(N_week*N_geo+1)]-theta0[(N_week*N_geo+1)]))
        # for (h in 1){
          if (runif(1,0,1) <= r){
            theta0[J] <- Ts[J]  # if accept, keep new parameter value
            L1 <- Lint          # if accept, keep new lieklihood
          }
        # }
      }else{
        for (h in 1:N_geo){
          if (runif(1,0,1) <= r[h]){
            theta0[J[h]] <- Ts[J[h]]  # if accept, keep new parameter value
            L1[h] <- Lint[h]          # if accept, keep new lieklihood
          }
        }
      }
      L[i,J] <- L1    # store final likelihood
    }
    
    thetas[i,] <- theta0  # store final parameter values for this iteration
  }
  
  ####
  res <- list(theta = thetas, logL = L)
  return(res)
}
