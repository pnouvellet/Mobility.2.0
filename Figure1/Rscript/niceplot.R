niceplot <- function(i,res_summary, res_base,inputs){
  
  fig.name=paste0('../../../../Dropbox (SPH Imperial College)/Mobility.2.0.Save/Saved_file/figure1/figs/',Mdata,'/',country[i],'.png')
  png(filename = fig.name, width = 9, height=6, units="in", res = 300)

  xlim <- c(as.Date('15-02-2020',format='%d-%m-%Y'),range(inputs$D$dates)[2])
  layout(matrix(c(1,1,3,3,3,2,2,3,3,3),2,5,byrow = TRUE))
  # plot mobility
  plot(inputs$D$dates,inputs$M_process_mat[[Mdata]][,k]/100,#main = country[i],
       xlim = xlim,
       ylim = c(0,max(c(1,inputs$M_process_mat[[Mdata]][,k]/100),na.rm=TRUE)),
       xlab = '',ylab = 'prop. mobility',pch=16,bty = 'n')
  lines(inputs$D$dates,inputs$M_process_mat[[Mdata]][,k]/100,col='darkorchid3',lwd=2)
  # lines(mob3[,i+1],col='red')
  abline(h=1,lwd=2,col='darkolivegreen4',lty=2)
  mtext('a)', side = 3, line = 1.5, adj = -.15, cex = 1)
  
  # plot Rt
  f <- 1:nrow(res_base[[si]]$resEpi$S_below.5) #which(res_base$resEpi$CV_below.2[,i+1])
  x <- res_base[[si]]$resEpi$median_R$dates[f] # (epi_res$t_end2 )[f]
  plot(x,res_base[[si]]$resEpi$median_R[f,k+1],
       xlim=xlim,
       ylim = c(0,5),
       type = 'l',col='black',lwd=2,bty = 'n',
       # main =country[i], 
       xlab = '',ylab = TeX('R_t'))
  mtext('b)', side = 3, line = 1.5, adj = -.15, cex = 1)
  
  polygon(c(x,rev(x)),
          c(res_base[[si]]$resEpi$low_R[f,k+1],rev(res_base[[si]]$resEpi$up_R[f,k+1])), 
          col = rgb(0,0,0,.2), border = FALSE )
  abline(h=(1),lwd=2,col='darkolivegreen4',lty=2)
  
  # plot new estimates
  f <- which(res_summary$results_full_Rt_D2$median_R[,i+1] != 0)
  lines(res_summary$results_full_Rt_D2$median_R$dates[f],
        res_summary$results_full_Rt_D2$median_R[f ,i+1],
        type='l',#ylim = c(0,10),
        # xlim=c(0,N_days),
        col='blue',lwd=2)
  
  polygon(c(res_summary$results_full_Rt_D2$median_R$dates[f],
            rev(res_summary$results_full_Rt_D2$median_R$dates[f])),
          c(res_summary$results_full_Rt_D2$low_R[f,i+1],
            rev(res_summary$results_full_Rt_D2$up_R[f,i+1])),
          col=rgb(0,0,1,.2), border = NA)
  
  lines(res_summary$results_full_Rt_daily$median_R$dates[f],
        res_summary$results_full_Rt_daily$median_R[f ,i+1],
        type='l',#ylim = c(0,10),
        # xlim=c(0,N_days),
        col='red',lwd=2)
  
  polygon(c(res_summary$results_full_Rt_daily$median_R$dates[f],
            rev(res_summary$results_full_Rt_daily$median_R$dates[f])),
          c(res_summary$results_full_Rt_daily$low_R[f,i+1],
            rev(res_summary$results_full_Rt_daily$up_R[f,i+1])),
          col=rgb(1,0,0,.2), border = NA)
  
  if (Change == 1){
    lines(inputs$D$dates[round(res_m$tchange[1,i])],5.5,type = 'p',pch=16,col='black')
    lines(inputs$D$dates[round(res_m$tchange[2:3,i])],rep(5.5,2),lwd=1.5)
  }
  if (Change ==0){
    legend('topleft',legend = c(TeX('R_t^{D}'),TeX('R_t'),TeX('R_t^{obs}')), 
           col = c(rgb(0,0,1),rgb(1,0,0),rgb(0,0,0)),lwd=2,bty='n',seg.len = .5)
  }else{
    legend('topleft',legend = c(TeX('R_t^{D}'),TeX('R_t'),TeX('R_t^{obs}'),TeX('t_{change}')), 
           col = c(rgb(0,0,1),rgb(1,0,0),rgb(0,0,0),rgb(0,0,0)),lwd=2,
           pch = c(rep(NA,3),16),bty='n',seg.len = .5)
  }
  
  
  
  # plot Rt vs mob
 
  
  f <- which(inputs$D[,k+1] > 0)
  if (Change == 0){
    plot(res_summary$results_full_Rt_assumed_mob$median_R$mobility,
         res_summary$results_full_Rt_assumed_mob$median_R[ ,i+1],
         type='l',#ylim = c(0,10),
         # xlim=c(0,N_days),
         col=rgb(.4,.5,.1),lwd=2, 
         xlab = 'Prop. reduction in movement',
         ylab = TeX('R_t'),yaxt="n",
         xlim = c(-1,1),bty = 'n',
         ylim = c((0.5),(6.5)))
    
    
    mtext('c)', side = 3, line = -2, adj = -.05, cex = 1)
    axis(side =2, at = (c(0,1,3,6)), labels = c(0,1,3,6))
    
    abline(h=(1),lwd=2,col='darkolivegreen4',lty=2)
    polygon(c(res_summary$results_full_Rt_assumed_mob$median_R$mobility,
              rev(res_summary$results_full_Rt_assumed_mob$median_R$mobility)),
            (c(res_summary$results_full_Rt_assumed_mob$low_R[,i+1],
               rev(res_summary$results_full_Rt_assumed_mob$up_R[,i+1]))),
            col=rgb(.4,.5,.1,.2),border=NA)
  }
  
  if (Change == 1){
    plot(res_summary$results_full_Rt_assumed_mob$median_R$mobility,
         res_summary$results_full_Rt_assumed_mob$median_R[ ,i+1],
         type='l',#ylim = c(0,10),
         # xlim=c(0,N_days),
         col=rgb(.8,.3,0),lwd=2, 
         xlab = 'Prop. reduction in movement',
         ylab = TeX('R_t'),yaxt="n",
         xlim = c(-1,1),bty = 'n',
         ylim = c((0.5),(6.5)))
    
    
    mtext('c)', side = 3, line = -2, adj = -.05, cex = 1)
    axis(side =2, at = (c(0,1,3,6)), labels = c(0,1,3,6))
    
    abline(h=(1),lwd=2,col='darkolivegreen4',lty=2)
    polygon(c(res_summary$results_full_Rt_assumed_mob$median_R$mobility,
              rev(res_summary$results_full_Rt_assumed_mob$median_R$mobility)),
            (c(res_summary$results_full_Rt_assumed_mob$low_R[,i+1],
               rev(res_summary$results_full_Rt_assumed_mob$up_R[,i+1]))),
            col=rgb(.8,.3,0,.2),border=NA)
    
    lines(res_summary$results_full_Rt_assumed_mob2$median_R$mobility,
         res_summary$results_full_Rt_assumed_mob2$median_R[ ,i+1],
         type='l',#ylim = c(0,10),
         # xlim=c(0,N_days),
         col=rgb(.4,.5,.1),lwd=2, 
         xlab = 'Prop. reduction in movement',
         ylab = TeX('R_t'),yaxt="n",
         xlim = c(0,1),bty = 'n',
         ylim = c((0.5),(6.5)))
    
    polygon(c(res_summary$results_full_Rt_assumed_mob2$median_R$mobility,
              rev(res_summary$results_full_Rt_assumed_mob2$median_R$mobility)),
            (c(res_summary$results_full_Rt_assumed_mob2$low_R[,i+1],
               rev(res_summary$results_full_Rt_assumed_mob2$up_R[,i+1]))),
            col=rgb(.4,.5,.1,.2),border=NA)
  }
  
  # f2 <- range(1-res_summary$results_meff$median_meff[f,i+1])
  # f3 <- which((res_summary$results_full_Rt_assumed_mob$median_R$mobility > f2[1]) & 
  #               (res_summary$results_full_Rt_assumed_mob$median_R$mobility < f2[2]))
  # lines(res_summary$results_full_Rt_assumed_mob$median_R$mobility[f3],
  #       res_summary$results_full_Rt_assumed_mob$median_R[f3,i+1],
  #       type='l',#ylim = c(0,10),
  #       # xlim=c(0,N_days),
  #       col=rgb(1,0,1),lwd=2)
  # 
  # polygon(c(res_summary$results_full_Rt_assumed_mob$median_R$mobility[f3],
  #           rev(res_summary$results_full_Rt_assumed_mob$median_R$mobility[f3])),
  #         (c(res_summary$results_full_Rt_assumed_mob$low_R[f3,i+1],
  #            rev(res_summary$results_full_Rt_assumed_mob$up_R[f3,i+1]))),
  #         col=rgb(1,0,1,.2),border=NA)
  
  #epiestim
  f <- seq(4,nrow(res_base[[si]]$resEpi$S_below.5),by=7) 
  which(res_base[[si]]$resEpi$S_below.5[,i+1])
  x <- 1-res_summary$results_meff$median_meff[f,i+1]-.5e-2
  y <- res_base[[si]]$resEpi$median_R[f,k+1]
  yplus <- res_base[[si]]$resEpi$up_R[f,k+1]
  yminus <- res_base[[si]]$resEpi$low_R[f,k+1]
  n_est <- length(x)
  
  if(Change == 0){
    # f_new <-  n_est-round(as.numeric(substr(Tb_best$t_change[i], 1,2))/7)
    errbar(x = x,
           y = y,
           yplus = yplus,
           yminus = yminus,
           col = c(rep(rgb(0,0,0),n_est-0),rep('red4',0)),
           errbar.col = c(rep(rgb(0,0,0),n_est-0),rep('red4',0)),
           add = TRUE)
    
    
    legend('topright',legend = c(TeX('R_t predicted (with mobility)'),
                                 # TeX('R_t within range of observed mobility'),
                                 TeX('past R_t^{obs}')),
           col = c(rgb(.4,.5,.1,.2),rgb(0,0,0)),lwd=2,bty='n')
  }else{
    f_new <-  n_est-round(as.numeric(substr(Tb_best$t_change[i], 1,3))/7)
    errbar(x = x,
           y = y,
           yplus = yplus,
           yminus = yminus,
           col = c(rep(rgb(0,0,0),n_est-f_new),rep('red4',f_new)),
           errbar.col = c(rep(rgb(0,0,0),n_est-f_new),rep('red4',f_new)),
           add = TRUE)
    
    
    legend('topright',legend = c(TeX('past R_t predicted (with mobility)'),
                                 TeX('new R_t predicted (with mobility)'),
                                 # TeX('R_t within range of observed mobility'),
                                 TeX('past R_t^{obs}'),
                                 TeX('past R_t^{obs} associated with new Th')),
           col = c(rgb(.8,.3,0,.2),rgb(.4,.5,.1,.2),'black','red4'),lwd=2,bty='n')
  }
  
  
  # if(Mdata == 'A'){
    mtext(paste0(country[i],' - ',Mdata), side = 3,  line = 1.5, outer = FALSE,adj = -.15,font=2)
  # }else{
    # mtext(paste0(country[i],' (Google)'), side = 3,  line = 1.5, outer = FALSE,adj = -.15,font=2)
    
  # }
  
  dev.off()
  
}