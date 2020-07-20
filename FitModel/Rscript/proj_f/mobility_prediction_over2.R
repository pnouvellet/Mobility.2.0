mobility_prediction_over2 <- function(n_pred,rep_sim){
  
  D_pred <- array(NA,dim = c(n_pred,N_geo,rep_sim))
  # new mobility
  
  H_pred <- matrix(0,nrow = N_days+n_pred, ncol = N_days+n_pred)
  for (i in 1:(N_days+n_pred)){
    f <- max(c(1,i-delta_id$SItrunc))
    H_pred[i,f:i] <- rev(delta_id$dist)[((delta_id$SItrunc+1)-(i-f)):(delta_id$SItrunc+1)]
    if (i>1) H_pred[i,f:i] <-  H_pred[i,f:i]/sum( H_pred[i,f:i])
  }
  
  
  # Rt_D2
  Rt_daily_pred <- array(NA,dim = c(nrow(M_pred),N_geo,rep_sim))
  Rt_D2_pred <- Rt_daily_pred
  if (rep_sim<=rep){
    k <- sample(x = 1:nrow(res$theta), size = rep_sim, replace = FALSE)  
  }else{
    k <- sample(x = 1:nrow(res$theta), size = rep_sim, replace = TRUE)  
  }
  for (j in 1:rep_sim){
    R_daily <- Rt_fun(theta = res$theta[k[j],], x = M_pred )
    Rt_daily_pred[,,j] <- R_daily
    Rt_D2_pred[,,j] <- H_pred %*% R_daily
    over <- res$theta[k[j],2*N_geo+1]
  }
  
  
  for (c in 1:N_geo){
    # project forward
    ws <- rev(SI$dist)
    
    D_temp <- matrix(c(mD[,c],rep(0,n_pred)),nrow(M_pred),rep_sim,byrow = FALSE)
    for (i in (N_days+(1:n_pred))){
      # lambda=t(D_temp[(i-SI$SItrunc):i,])%*%ws
      f <- max(c(1,(i-SI$SItrunc)))
      lambda=t(D_temp[f:i,])%*%ws[((SI$SItrunc+1)-(i-f)):(SI$SItrunc+1)]
      # D_temp[i,]=rpois(rep_sim,Rt_D2_pred[i,c,]*lambda)
      lambda[lambda>1e5] <- 1e5
      lambda[lambda==0] <- 1e-20
      D_temp[i,]=rnbinom(n = rep_sim, mu = Rt_D2_pred[i,c,]*lambda, 
                         size = sqrt((Rt_D2_pred[i,c,]*lambda)) * over) #over)  #(Rt_D2_pred[i,c,]*lambda)/over)
      
    }
    D_pred[,c,] <- D_temp[(N_days+(1:n_pred)),]
  }
  
  res <- list(D_pred = D_pred,
              Rt_D2_pred = Rt_D2_pred[(N_days+(1:n_pred)),,],
              Rt_daily_pred = Rt_daily_pred[(N_days+(1:n_pred)),,])
}