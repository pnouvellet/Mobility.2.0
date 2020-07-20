process_mobility <- function(d_m){
  
  dates_c <- as.character(substr(names(d_m[,-c(1:3)]),start = 2, stop = 11))
  dates <- as.Date(dates_c,format='%Y.%m.%d')
  
  n_day_m <- round(length(dates)/7)*7
  # dates <- dates[1:n_day_m]
  
  mob_raw <- data.frame(dates = dates,
                        matrix(NA,n_day_m,N_geo))
  names(mob_raw) <- c('dates',country)
  
  mob_combined_smooth <- data.frame(dates = dates,
                                    matrix(NA,n_day_m,N_geo))
  names(mob_combined_smooth) <- c('dates',country)
  # mob3 <- data.frame(dates = dates,
  #                   matrix(NA,n_day_m,N_geo))
  
  for (i  in 1:N_geo){
    f <- which(d_m$region %in% country[i])
    m <- d_m[f,]
    x <- colSums(m[,-c(1:3)])
    mob_raw[,i+1] <- x[1:n_day_m]
    
    x2 <- rowSums(matrix(mob_raw[,i+1],n_day_m/7,7,by = TRUE))/7
    mob_combined_smooth[4+seq(1,(n_day_m/7)*7,by=7),i+1] <- x2
    mob_combined_smooth[,i+1] <- na.approx(mob_combined_smooth[,i+1],rule = 2)
    
    # rescaling
    max_mob <- max(mob_combined_smooth[,i+1])
    mob_combined_smooth[,i+1] <- mob_combined_smooth[,i+1]/max_mob
    mob_raw[,i+1] <- mob_raw[,i+1]/max_mob
  }
  
  return(list(mob_raw = mob_raw,
              mob_combined_smooth = mob_combined_smooth ))
  
}