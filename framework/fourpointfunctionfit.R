require(parallel)
require(boot)
require(hadron)
source("framework/utils.R")
source("framework/parameters.R")

load("obskaon.RData")
load("obspion.RData")

pionkaonmass_summary <- data.frame()
obslistpionkaon <- list()
for (beta in betas){
  betastr <- make.betastr(beta)
  obslistpionkaon[[betastr]] <- list()
  for (mul in muls[[betastr]]){
    mulstr <- make.mulstr(mul)
    obslistpionkaon[[betastr]][[mulstr]] <-list()
    for (L in Ls){
      Lstr <- make.Lstr(L)
      obslistpionkaon[[betastr]][[mulstr]][[Lstr]] <-list()
      for (mus in  mussindex[[betastr]]){  
        musstr <- make.musstr(mus/100000)
   #     print(as.integer(as.double(mus*100000)))
        new<-try(readtextcf(file=sprintf("/home/pittler/Bonn_postdoc/pik_scattering/pik_charged_A1_TP0_00_%s%03d_%d_amus_%4d.dat",latticelist[[match(beta,betas)]],mul*10000,L,mus), T=2*L, sym=TRUE, path="", skip=1, check.t=0))
        if( any(class(new) == "try-error" ) ) {
          next
        }
 
        tikzfiles <- tikz.init(sprintf("quick_analysis_pik%s%03d_%d_%4d",latticelist[[match(beta,betas)]],mul*10000,L,mus), width=4.5, height=4.5)

 
        tmin <- fitrangepionkaon[[betastr]][[mulstr]][[Lstr]][2]
        tmax <- fitrangepionkaon[[betastr]][[mulstr]][[Lstr]][3]
        lmin <- fitrangepionkaon[[betastr]][[mulstr]][[Lstr]][1]
        lmax <- L


        obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]] <- list()
        obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]] <- list()
        obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]] <- list()

        for (useCov in c(FALSE,TRUE)) {
          useCovstr <- make.covstr(useCov)

          mkaon <- obslistkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][[useCovstr]][["M_ka_corr"]]$val
          mpion <- obslistpion[[betastr]][[mulstr]][[Lstr]][[useCovstr]][["M_pi_corr"]]$val
          edeltaE <- exp((mkaon - mpion)*(rep(1,L+1)))
          pollution <- exp((mkaon - mpion)*(seq(0,L)))*exp(-mkaon*2*L) + exp(-(mkaon - mpion)*(seq(0,L)))*exp(-mpion*2*L)
     
          newb<- bootstrap.cf(cf=new,boot.R=boot.R,boot.l=blocklength[[betastr]][[mulstr]][[Lstr]],seed=551,endcorr=TRUE)
          weightedcorr_E1 <- edeltaE * newb$cf.tsboot$t0
          weightedcorr_E2 <- pollution[seq(1,L)]/pollution[seq(2,L+1)]*newb$cf.tsboot$t0[seq(2,L+1)]

          correlation_E1 <- newb$cf.tsboot$t0[seq(1,L)]-weightedcorr_E1[seq(2,L+1)] 
          correlation_E2 <- newb$cf.tsboot$t0[seq(1,L)]-weightedcorr_E2[seq(1,L)]

          mkaon_bs <- obslistkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["M_ka_corr"]]$fit_bs
          mpion_bs <- obslist[[betastr]][[mulstr]][[Lstr]][["M_pi_corr"]]$fit_bs

          cd <- lapply(X=1:boot.R,FUN=function(rw_idx){newb$cf.tsboot$t[rw_idx,seq(1,L)]-exp((mkaon_bs[rw_idx,1] - mpion_bs[rw_idx,1])*(rep(1,L)))*newb$cf.tsboot$t[rw_idx,seq(2,L+1)];})
          fit_bs_E1 <- do.call(rbind, cd)

          cd <- lapply(X=1:boot.R,FUN=function(rw_idx){exp((mkaon_bs[rw_idx,1] - mpion_bs[rw_idx,1])*(seq(0,L)))*exp(-mkaon_bs[rw_idx,1]*2*L)+exp(-(mkaon_bs[rw_idx,1] - mpion_bs[rw_idx,1])*(seq(0,L)))*exp(-mpion_bs[rw_idx,1]*2*L)})
        fit_bs_pollution <- do.call(rbind, cd)

          cd <- lapply(X=1:boot.R,FUN=function(rw_idx){newb$cf.tsboot$t[rw_idx,seq(1,L)]-fit_bs_pollution[rw_idx,seq(1,L)]/fit_bs_pollution[rw_idx,seq(2,L+1)]*newb$cf.tsboot$t[rw_idx,seq(2,L+1)];})
          fit_bs_E2 <- do.call(rbind, cd)

    
          pvaluefit_E1 <- vector()
          fitresultsrange_E1 <- vector()
          fitresultsrange1_E1 <- vector()
          fitresultsrange2_E1 <- vector()
          fitresultsbs_E1 <- vector()
          fitresultsbs1_E1 <- vector()
          fitresultsbs2_E1 <- vector()
          fitresultserr_E1 <- vector()
          fitresultschisq_E1 <- vector ()

          pvaluefit_E2 <- vector()
          fitresultsrange_E2 <- vector()
          fitresultsrange1_E2 <- vector()
          fitresultsbs_E2 <- vector()
          fitresultsbs1_E2 <- vector()
          fitresultserr_E2 <- vector()
          fitresultschisq_E2 <- vector ()


          for (t1 in tmin:(tmax-lmin)){
            if (useCov == TRUE){
              t2min <- t1+lmin
              t2max <- tmax
            }
            else{
              t2min <- tmax
              t2max <- tmax
            }
            for (t2 in t2min:t2max){
              if (useCov == TRUE){
                 Cinv<-(solve(cov(fit_bs_E1[,t1:t2])))
              } else if (useCov == FALSE) {
                 cd<-lapply(X=t1:t2,FUN=function(rw_idx){sd(fit_bs_E1[,rw_idx]);})
                 fit_bs <- do.call(rbind, cd)
                 Cinv<-diag(as.vector(1/fit_bs^2))
              }
            
      #Cinv<-diag(1/sd(newb$cf.tsboot$t[,17:31])^2)
              pars.init<-initpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]]
              mfitglobal<- optim(par=pars.init,
                 method="BFGS",
                 # lambda function  
                 fn=chisqfn_4pt,
                 # here one could pass the first derivative of the chisqfn
                 # to make the fit behave better
                 gr=NULL,#chisqfnderivative,
                 # get some debug output from the algorithm
                 control=list(trace=TRUE,maxit=1000000),
                 # finally, the remaining arguments for the chisqfn
                 t=seq(t1-1,t2-1),
                 T=2*L,
                 corrdata =correlation_E1[t1:t2],
                 fitfn=fourpointfit_E1,
                 mpion=mpion,
                 mkaon=mkaon,
                 Cinv=Cinv
                )
              print(mfitglobal$par)
              mfitglobalfinal<- optim(par=mfitglobal$par,
                 method="BFGS",
                 # lambda function  
                 fn=chisqfn_4pt,
                 # here one could pass the first derivative of the chisqfn
                 # to make the fit behave better
                 gr=NULL,#chisqfnderivative,
                 # get some debug output from the algorithm
                 control=list(trace=TRUE,maxit=1000000),
                 # finally, the remaining arguments for the chisqfn
                 t=seq(t1-1,t2-1),
                 T=2*L,
                 corrdata =correlation_E1[t1:t2],
                 fitfn=fourpointfit_E1,
                 mpion=mpion,
                 mkaon=mkaon,
                 Cinv=Cinv
                )

              if (useCov == TRUE){
                pvaluefit_E1 <- rbind(pvaluefit_E1,pchisq(mfitglobal$value, df=t2-t1-2, lower.tail=FALSE))
              }
              fitresultsrange_E1 <- rbind(fitresultsrange_E1,mfitglobal$par[2])
              fitresultsrange1_E1 <- rbind( fitresultsrange1_E1, mfitglobal$par[1])
              fitresultsrange2_E1 <- rbind( fitresultsrange2_E1, mfitglobal$par[3])
              fitresultschisq_E1 <- rbind(fitresultschisq_E1, mfitglobal$value)
            
              res_lst <- mclapply(X = 1:boot.R,
                FUN= function(rw_idx){
                  mfit<- optim(par=mfitglobal$par,
                   method="BFGS",
                   # lambda function  
                   fn=chisqfn_4pt,
                   # here one could pass the first derivative of the chisqfn
                   # to make the fit behave better
                   gr=NULL,#chisqfnderivative,
                   # get some debug output from the algorithm
                   control=list(trace=FALSE,maxit=1000000),
                   # finally, the remaining arguments for the chisqfn
                   t=seq(t1-1,t2-1),
                   T=2*L,
                   corrdata = fit_bs_E1[rw_idx,t1:t2],
                   fitfn=fourpointfit_E1,
                   mpion=mpion_bs[rw_idx],
                   mkaon=mkaon_bs[rw_idx],
                   Cinv=Cinv
                  )
                  return(mfit$par)
                }
              )
              fit_bs <- do.call(rbind, res_lst)
              
              fitresultsbs1_E1 <- cbind(fitresultsbs1_E1, fit_bs[,1])
              fitresultsbs2_E1 <- cbind(fitresultsbs2_E1, fit_bs[,3])
              fitresultsbs_E1 <- cbind(fitresultsbs_E1, fit_bs[,2])
              fitresultserr_E1 <- rbind(fitresultserr_E1, sd(fit_bs[,2]))  
            

              if (useCov == TRUE){
                Cinv<-(solve(cov(fit_bs_E2[,t1:t2])))
              } else if (useCov == FALSE) {
                cd<-lapply(X=t1:t2,FUN=function(rw_idx){sd(fit_bs_E2[,rw_idx]);})
                fit_bs <- do.call(rbind, cd)
                Cinv<-diag(as.vector(1/fit_bs^2))
              }

      #Cinv<-diag(1/sd(newb$cf.tsboot$t[,17:31])^2)
              pars.init<-c(initpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][1],initpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][2])
              mfitglobal<- optim(par=pars.init,
                 method="BFGS",
                 # lambda function  
                 fn=chisqfn_4pt,
                 # here one could pass the first derivative of the chisqfn
                 # to make the fit behave better
                 gr=NULL,#chisqfnderivative,
                 # get some debug output from the algorithm
                 control=list(trace=TRUE,maxit=1000000),
                 # finally, the remaining arguments for the chisqfn
                 t=seq(t1-1,t2-1),
                 T=2*L,
                 corrdata =correlation_E2[t1:t2],
                 fitfn=fourpointfit_E2,
                 mpion=mpion,
                 mkaon=mkaon,
                 Cinv=Cinv
                )

              print(mfitglobal$par)
              if (useCov == TRUE){
                pvaluefit_E2 <- rbind(pvaluefit_E2,pchisq(mfitglobal$value, df=t2-t1-1, lower.tail=FALSE))
              }
              fitresultsrange_E2 <- rbind(fitresultsrange_E2,mfitglobal$par[2])
              fitresultsrange1_E2 <- rbind( fitresultsrange1_E2, mfitglobal$par[1])
              fitresultschisq_E2 <- rbind(fitresultschisq_E2, mfitglobal$value)

              res_lst <- mclapply(X = 1:boot.R,
                FUN= function(rw_idx){
                  mfit<- optim(par=pars.init,
                     method="BFGS",
                     # lambda function  
                     fn=chisqfn_4pt,
                     # here one could pass the first derivative of the chisqfn
                     # to make the fit behave better
                     gr=NULL,#chisqfnderivative,
                     # get some debug output from the algorithm
                     control=list(trace=FALSE,maxit=1000000),
                     # finally, the remaining arguments for the chisqfn
                     t=seq(t1-1,t2-1),
                     T=2*L,
                     corrdata = fit_bs_E2[rw_idx,t1:t2],
                     fitfn=fourpointfit_E2,
                     mpion=mpion_bs[rw_idx],
                     mkaon=mkaon_bs[rw_idx],
                     Cinv=Cinv
                  )
                  return(mfit$par)
                }
              )
              fit_bs <- do.call(rbind, res_lst)
              fitresultsbs1_E2 <- cbind(fitresultsbs1_E2, fit_bs[,1])
              fitresultsbs_E2 <- cbind(fitresultsbs_E2, fit_bs[,2])
              fitresultserr_E2 <- rbind(fitresultserr_E2, sd(fit_bs[,2]))
           
            }#t2
          }#t1
          if (useCov == TRUE){
            minstat <- min(fitresultserr_E1)
            weightmeff<- ((1-2*(pvaluefit_E1-0.5)^2)*minstat/fitresultserr_E1)^2
            finalres <- weighted.median(x=fitresultsrange_E1,w=weightmeff)
            cd <- lapply(X=1:boot.R,FUN=function(rw_idx){weighted.median(x=fitresultsbs_E1[rw_idx,],w=weightmeff)})
            fit_bs <- do.call(rbind, cd)
            normalizedweight_E1 <- (weightmeff/sum(weightmeff)) 
          }
          if (useCov == FALSE){
            x <- 1
            listdev <- vector()

            for (t1 in tmin:(tmax-lmin)){
              for (t2 in tmax){
                listt <- rep(1,t2-t1+1)
                minstat <- min(fitresultserr)
                listdev <- rbind(listdev,sum(abs(correlation_E1[t1:t2]/fourpointfit_E1((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E1[x],fitresultsrange_E1[x],fitresultsrange2_E1[x]))-listt))/length(seq(t1-1,t2-1)))
#              listdev <- rbind(listdev,abs(sum(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x])))/length(seq(t1-1,t2-1))-1))
#              listdev <- rbind(listdev,temporary)
                x <- x + 1
              }
            }
            print(listdev)
            finalres <- fitresultsrange_E1[which.min(listdev)]
            fit_bs <- fitresultsbs_E1[,which.min(listdev)]

          }

          obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]] <-list(val=finalres,dval=sd(fit_bs),fit_bs=fit_bs)
           
          if (useCov == TRUE ){
            minstat <- min(fitresultserr_E2)
            weightmeff<- ((1-2*(pvaluefit_E2-0.5)^2)*minstat/fitresultserr_E2)^2
            finalres <- weighted.median(x=fitresultsrange_E2,w=weightmeff)
            cd <- lapply(X=1:boot.R,FUN=function(rw_idx){weighted.median(x=fitresultsbs_E2[rw_idx,],w=weightmeff)})
            fit_bs <- do.call(rbind, cd)
            normalizedweight_E2 <- (weightmeff/sum(weightmeff))
          }
          if (useCov == FALSE){
            x <- 1
            listdev <- vector()

            for (t1 in tmin:(tmax-lmin)){
              for (t2 in tmax){
                listt <- rep(1,t2-t1+1)
                minstat <- min(fitresultserr)
                listdev <- rbind(listdev,sum(abs(correlation_E2[t1:t2]/fourpointfit_E2((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E2[xdd],fitresultsrange_E2[xdd]))-listt))/length(seq(t1-1,t2-1)))
#              listdev <- rbind(listdev,abs(sum(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x])))/length(seq(t1-1,t2-1))-1))
#              listdev <- rbind(listdev,temporary)
                x <- x + 1
              }
            }
            print(listdev)
            finalres <- fitresultsrange_E2[which.min(listdev)]
            fit_bs <- fitresultsbs_E2[,which.min(listdev)]

          }

          obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][[useCovstr]] <-list(val=finalres,dval=sd(fit_bs),fit_bs=fit_bs)

          pionkaonmass_summary <- rbind( pionkaonmass_summary,
                                    data.frame(
                                           beta=betastr,
                                           mul=mulstr,
                                           L=Lstr,
                                           mus=musstr,
                                           useCov=useCovstr,
                                           valE1=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val,
                                           dvalE1=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$dval,
                                           valE2=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][[useCovstr]]$val,
                                           dvalE2=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][[useCovstr]]$dval
                                    )
                                  )

#Compute the effective mass with E2

          fn <- function(m, t, T, corr_E2, mkaon,mpion) {
             pollutiont   <- exp(-mpion *T)*exp(-(mkaon-mpion)*(t-1))+exp(-mkaon *T)*exp((mkaon-mpion)*(t-1))
             pollutiontp1 <- exp(-mpion *T)*exp(-(mkaon-mpion)*(t+0))+exp(-mkaon *T)*exp((mkaon-mpion)*(t+0))
             pollutiontp2 <- exp(-mpion *T)*exp(-(mkaon-mpion)*(t+1))+exp(-mkaon *T)*exp((mkaon-mpion)*(t+1))
             return(corr_E2[(t+1)]/corr_E2[t] - (exp(-m * (t + 0)) + exp(-m *(T - (t + 0))) - pollutiontp1/pollutiontp2 * (exp(-m * (t + 1 )) + exp(-m * (T - (t + 1)))))/(exp(-m * (t -1)) + exp(-m *(T - (t - 1))) - pollutiont/pollutiontp1 * (exp(-m * (t + 0 )) + exp(-m * (T - (t + 0))))))
          }

          cd<-lapply(X=tmin:(tmax-1),FUN=function(rw_idx) {uniroot(fn, c(0,1),extendInt="yes", trace=1,t=rw_idx,T=2*L,corr_E2=correlation_E2,mkaon=mkaon,mpion=mpion)$root})
          effE2 <- do.call(rbind,cd)
          print(effE2)

          cd4 <-lapply(X=1:boot.R, FUN=function(rw_idx2){cd2<-lapply(X=tmin:(tmax-1),FUN=function(rw_idx) {uniroot(fn, c(0,1),extendInt="yes", trace=1,t=rw_idx,T=2*L,corr_E2=fit_bs_E2[rw_idx2,],mkaon=mkaon_bs[rw_idx2,],mpion=mpion_bs[rw_idx2,])$root});fitss<-do.call(cbind,cd2)})

          cd8<-lapply(X=1:(tmax-tmin),FUN=function(rw_idx2){cd6<-lapply(X=1:boot.R,FUN=function(rw_idx){cd4[[rw_idx]][rw_idx2]});fitss<-do.call(rbind,cd6);sd(fitss)})
          effE2_error<-do.call(rbind,cd8)
          print(effE2_error)

          xdd <- 1

          for (t1 in tmin:(tmax-lmin)){

            if (useCov == TRUE){
              t2min <- t1+lmin
              t2max <- tmax
            }
            else{
              t2min <- tmax
              t2max <- tmax
            }

            for (t2 in t2min:t2max){
              plotwitherror(
                x=(tmin-1):(tmax-2),
                y=effE2,
                dy=effE2_error,
                ylab="$a E_{\\pi K}^\\mathrm{eff}(t)$",
                xlab="$t/a$",
                main=sprintf("$\\beta %1.2f~a \\mu_l = %.4f~a \\mu_s = %d\\newline L=%d~C=%d~t_1=%d~t_2=%d$", beta,as.numeric(mul),as.numeric(mus),L,useCov,t1,t2),
                ylim=c(effE2[1]-5*effE2_error[1],effE2[1]+5*effE2_error[1]))
              legend(x="bottomleft",
                bty='n',
                pch=NA,
                lty=1,
                col=c(rep("black",2),rep("blue",2),rep("red",1),rep("green",1)),
                legend=c(sprintf("$a E_{\\pi K}^\\mathrm{E1} = %s $",
                            tex.catwitherror(x=fitresultsrange_E1[xdd],
                                             dx=fitresultserr_E1[xdd],
                                             digits=3,
                                             with.dollar=FALSE)),
                         sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                            fitresultschisq_E1[xdd]/(t2-t1+1-3)),
                         sprintf("$a E_{\\pi K}^\\mathrm{E2} = %s $",
                            tex.catwitherror(x=fitresultsrange_E2[xdd],
                                             dx=fitresultserr_E2[xdd],
                                             digits=3,
                                             with.dollar=FALSE)),
                         sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                            fitresultschisq_E2[xdd]/(t2-t1+1-2)),
                         sprintf("$a E_{\\pi K}^\\mathrm{weighted~median,E1} = %s $",
                            tex.catwitherror(x=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val,
                                             dx=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$dval,
                                             digits=3,
                                             with.dollar=FALSE)),
                         sprintf("$a E_{\\pi K}^\\mathrm{weighted~median,E2} = %s $",
                            tex.catwitherror(x=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][[useCovstr]]$val,
                                             dx=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][[useCovstr]]$dval,
                                             digits=3,
                                             with.dollar=FALSE))
                    )
              )
              pcol <- col2rgb("blue", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              rect(xleft=t1-1,
                   xright=t2-1,
                   ybottom=fitresultsrange_E1[xdd]-fitresultserr_E1[xdd],
                   ytop=fitresultsrange_E1[xdd]+fitresultserr_E1[xdd],
                   col=pcol,
                   border=NA
              )
              lines(x=c(t1,t2),
                 y=c(fitresultsrange_E1[xdd],fitresultsrange_E1[xdd]),
                 col="blue"
              )
              pcol <- col2rgb("red", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              rect(xleft=tmin,
                   xright=tmax,
                   ybottom=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val-obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$dval,
                   ytop=obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val+obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$dval,
                   col=pcol,
                   border=NA
              )
              lines(x=c(t1,t2),
                   y=c(obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val,obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][[useCovstr]]$val),
                   col="red"
              )

                 
              err<-NULL
              for (t in tmin:tmax){ temp<-lapply(X=1:boot.R,FUN=function(rw_idx){fit_bs_E1[rw_idx,t]/fourpointfit_E1(t-1,2*L,mpion_bs[rw_idx],mkaon_bs[rw_idx],c(fitresultsbs1_E1[rw_idx,xdd],fitresultsbs_E1[rw_idx,xdd],fitresultsbs2_E1[rw_idx,xdd]))});temp2 <-do.call(rbind, temp);err<-rbind(err,sd(as.vector(temp2)))}


              plotwitherror(x = (tmin-1):(tmax-1),
                            y = correlation_E1[tmin:tmax]/fourpointfit_E1((tmin-1):(tmax-1),2*L,mpion,mkaon,c(fitresultsrange1_E1[xdd],fitresultsrange_E1[xdd],fitresultsrange2_E1[xdd])),
                            dy = err,
                            ylab = "$C(t)/f(t)$",
                            xlab="$t/a$",
                            main=sprintf("$E1~weight=%4f~t_1=%d~t_2=%d$", normalizedweight_E1[xdd],t1,t2)
                           )
              pcol <- col2rgb("red", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              erma<-max(err)
              rect(xleft=t1-1,
                   xright=t2-1,
                   ybottom=min(correlation_E1[t1:t2]/fourpointfit_E1((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E1[xdd],fitresultsrange_E1[xdd],fitresultsrange2_E1[xdd]))-erma),
                   ytop=max(correlation_E1[t1:t2]/fourpointfit_E1((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E1[xdd],fitresultsrange_E1[xdd],fitresultsrange2_E1[xdd]))+erma),
                   col=pcol,
                   border=NA
              )

              err<-NULL
              for (t in tmin:tmax){ temp<-lapply(X=1:boot.R,FUN=function(rw_idx){fit_bs_E2[rw_idx,t]/fourpointfit_E2(t-1,2*L,mpion_bs[rw_idx],mkaon_bs[rw_idx],c(fitresultsbs1_E2[rw_idx,xdd],fitresultsbs_E2[rw_idx,xdd]))});temp2 <-do.call(rbind, temp);err<-rbind(err,sd(as.vector(temp2)))}


              plotwitherror(x = (tmin-1):(tmax-1),
                            y = correlation_E2[tmin:tmax]/fourpointfit_E2((tmin-1):(tmax-1),2*L,mpion,mkaon,c(fitresultsrange1_E2[xdd],fitresultsrange_E2[xdd])),
                            dy = err,
                            ylab = "$C(t)/f(t)$",
                            xlab="$t/a$",
                            main=sprintf("$E2~weight=%4f~t_1=%d~t_2=%d$", normalizedweight_E2[xdd],t1,t2)
                       )
              pcol <- col2rgb("red", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              erma <- max(err)
              rect(xleft=t1-1,
                   xright=t2-1,
                   ybottom=min(correlation_E2[t1:t2]/fourpointfit_E2((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E2[xdd],fitresultsrange_E2[xdd]))-erma),
                   ytop=max(correlation_E2[t1:t2]/fourpointfit_E2((t1-1):(t2-1),2*L,mpion,mkaon,c(fitresultsrange1_E2[xdd],fitresultsrange_E2[xdd]))+erma),
                   col=pcol,
                   border=NA
              )



              xdd <- xdd + 1

            }#t2
          }#t1
        }#useCovstr
        tikz.finalize(tikzfiles)
      }#mus
    }#L
  }#mul
}#beta
save(obslistpionkaon,
     file = "obspionkaon.RData")
save(pionkaonmass_summary,
     file = "obspionkaon_summary.RData")
