require(parallel)
require(boot)
require(hadron)
source("framework/utils.R")
source("framework/parameters.R")

kaonmass_summary<-data.frame()
obslistkaon <- list()

for (beta in betas){
  betastr <- make.betastr(beta)
  obslistkaon[[betastr]] <- list()
  for (mul in muls[[betastr]]){
    mulstr <- make.mulstr(mul)
    obslistkaon[[betastr]][[mulstr]] <-list()
    for (L in Ls){
      Lstr <- make.Lstr(L)
      obslistkaon[[betastr]][[mulstr]][[Lstr]] <- list()
      for (mus in mussindex[[betastr]] ){  
        Lstr <- make.Lstr(L)
        musstr <- make.musstr(mus/100000)
   #     print(as.integer(as.double(mus*100000)))
        new<-try(readtextcf(file=sprintf("/home/pittler/Bonn_postdoc/pik_scattering/k_charged_p0_%s%03d_%d_amus_%4d.dat",latticelist[[match(beta,betas)]],mul*10000,L,mus), T=2*L, sym=TRUE, path="", skip=1, check.t=0, symmetrise=TRUE))
        if( any(class(new) == "try-error" ) ) {
          next
        }
        tikzfiles <- tikz.init(sprintf("quick_analysiskaon%s%03d_%d_%4d_abs",latticelist[[match(beta,betas)]],mul*10000,L,mus), width=4.5, height=4.5)

        obslistkaon[[betastr]][[mulstr]][[Lstr]][[musstr]] <- list()


 
        tmin <- fitrangekaon[[betastr]][[mulstr]][[Lstr]][2]
        tmax <- fitrangekaon[[betastr]][[mulstr]][[Lstr]][3]
        lmin <- fitrangekaon[[betastr]][[mulstr]][[Lstr]][1]
        lmax <- L

        newb<- bootstrap.cf(cf=new,boot.R=boot.R,boot.l=blocklength[[betastr]][[mulstr]][[Lstr]],seed=551,endcorr=TRUE)

        for (useCov in c( FALSE, TRUE )){

          pvaluefit <- vector()
          fitresultsrange <- vector()
          fitresultsrange1 <- vector()
          fitresultsbs <- vector()
          fitresultsbs1 <- vector()
          fitresultserr <- vector()
          fitresultschisq <- vector ()

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
                Cinv<-(solve(cov(newb$cf.tsboot$t[,t1:t2])))
              } else if (useCov == FALSE) {
                cd<-lapply(X=t1:t2,FUN=function(rw_idx){sd(newb$cf.tsboot$t[,rw_idx]);})
                fit_bs <- do.call(rbind, cd)
                Cinv<-diag(as.vector(1/fit_bs^2))
              }
  
              pars.init<-initkaon[[betastr]][[mulstr]][[Lstr]][[musstr]]
              mfitglobal<- optim(par=pars.init,
                 method="BFGS",
                 # lambda function  
                 fn=chisqfn,
                 # here one could pass the first derivative of the chisqfn
                 # to make the fit behave better
                 gr=NULL,#chisqfnderivative,
                 # get some debug output from the algorithm
                 control=list(trace=TRUE,maxit=1000000),
                 # finally, the remaining arguments for the chisqfn
                 t=seq(t1-1,t2-1),
                 T=2*L,
                 corrdata =newb$cf.tsboot$t0[t1:t2],
                 fitfn=coshfit,
                 Cinv=Cinv
              )
              #Doing a second fit to test the convergence
              #and use the results of the first as an input
              mfitglobalfinal<-  optim(par=mfitglobal$par,
                 method="BFGS",
                 # lambda function  
                 fn=chisqfn,
                 # here one could pass the first derivative of the chisqfn
                 # to make the fit behave better
                 gr=NULL,#chisqfnderivative,
                 # get some debug output from the algorithm
                 control=list(trace=TRUE,maxit=1000000),
                 # finally, the remaining arguments for the chisqfn
                 t=seq(t1-1,t2-1),
                 T=2*L,
                 corrdata =newb$cf.tsboot$t0[t1:t2],
                 fitfn=coshfit,
                 Cinv=Cinv
              )

              if (useCov == TRUE){
                pvaluefit <- rbind(pvaluefit,pchisq(mfitglobalfinal$value, df=t2-t1-1, lower.tail=FALSE))
              }

              fitresultsrange <- rbind(fitresultsrange,mfitglobalfinal$par[2])
              fitresultsrange1 <- rbind( fitresultsrange1, mfitglobalfinal$par[1])
              fitresultschisq <- rbind(fitresultschisq, mfitglobalfinal$value)
              
              res_lst <- mclapply(X = 1:boot.R,
                FUN= function(rw_idx){
                  mfit<- optim(par=mfitglobal$par,
                     method="BFGS",
                     # lambda function  
                     fn=chisqfn,
                     # here one could pass the first derivative of the chisqfn
                     # to make the fit behave better
                     gr=NULL,#chisqfnderivative,
                     # get some debug output from the algorithm
                     control=list(trace=FALSE,maxit=1000000),
                     # finally, the remaining arguments for the chisqfn
                     t=seq(t1-1,t2-1),
                     T=2*L,
                     corrdata =newb$cf.tsboot$t[rw_idx,t1:t2],
                     fitfn=coshfit,
                       Cinv=Cinv
                    )
                  return(mfit$par)
                  }
                )
                fit_bs <- do.call(rbind, res_lst)
                fitresultsbs1 <- cbind(fitresultsbs1, fit_bs[,1])
                fitresultsbs <- cbind(fitresultsbs, fit_bs[,2])
                fitresultserr <- rbind(fitresultserr, sd(fit_bs[,2]))  
              }#t2
            }#t1

          if( useCov == TRUE){
            minstat <- min(fitresultserr)
            weightlist <- 1:length(fitresultserr)
#Calculate the weight of the fit using the pvalue and the statistical error
#We calculate the weight factors only for the effective mass and apply it to get
#the weighted median
            weightmeff <- ((1-2*(pvaluefit-0.5)^2)*minstat/fitresultserr)^2
            finalres <- weighted.median(x=fitresultsrange,w=weightmeff)
            cd <- lapply(X=1:boot.R,FUN=function(rw_idx){weighted.median(x=fitresultsbs[rw_idx,],w=weightmeff)})
            fit_bs <- do.call(rbind, cd)
            normalizedweight <- (weightmeff/sum(weightmeff))
          }#useCov true
          if (useCov == FALSE){
            x <- 1
            listdev <- vector()
            for (t1 in tmin:(tmax-lmin)){
              for (t2 in tmax){
                listt <- rep(1,t2-t1+1)
                listdev <- rbind(listdev,sum(abs(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))-listt))/length(seq(t1-1,t2-1)))
                x <- x + 1
              }
            }
            print(listdev)
            finalres <- fitresultsrange[which.min(listdev)]
            fit_bs <- fitresultsbs[,which.min(listdev)]
          }

          useCovstr <- make.covstr(useCov)

          obslistkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][[useCovstr]] <- list()
          obslistkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][[useCovstr]][["M_ka_corr"]] <-list(val=finalres,dval=sd(fit_bs),fit_bs=fit_bs)


          kaonmass_summary <- rbind( kaonmass_summary,
                                    data.frame(
                                           beta=betastr,
                                           mul=mulstr,
                                           L=Lstr,
                                           mus=musstr,
                                           useCov=useCovstr,
                                           val=finalres,
                                           dval=sd(fit_bs)
                               )
                              )


          x <- 1
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
              listt <- rep(1,t2-t1+1)
              if (t2 == tmax){
                efmassfit <- fit.effectivemass(
                     cf=bootstrap.effectivemass(cf=newb,
                                          type='solve'),
                                          t1=t1-1,
                                          t2=t2-2,
                                          useCov=useCov
                 )
              }
              else{
                efmassfit <- fit.effectivemass(
                     cf=bootstrap.effectivemass(cf=newb,
                                          type='solve'),
                                          t1=t1-1,
                                          t2=t2-1,
                                          useCov=useCov
                 )
              }

              m_eff <- efmassfit$effmassfit

              plot(efmassfit,
                ylab="$a M_K^\\mathrm{eff}(t)$",
                xlab="$t/a$",
                main=sprintf("$\\beta %1.2f~a \\mu_l = %.4f~a \\mu_s = %d\\newline L=%d~Cov=%d~t_1=%d~t_2=%d$", beta,as.numeric(mul),as.numeric(mus),L,useCov,t1-1,t2-1),
                ylim=c(m_eff$t0[1]-5*m_eff$se,m_eff$t0[1]+5*m_eff$se))
                legend(x="bottomleft",
                bty='n',
                pch=NA,
                lty=1,
                col=c(rep("black",2),rep("blue",2),rep("red",1)),
                legend=c(sprintf("$a M_K^\\mathrm{eff} = %s $",
                             tex.catwitherror(x=m_eff$t0[1],
                                              dx=m_eff$se,
                                              digits=3,
                                              with.dollar=FALSE)),
                         sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                                              m_eff$t0[2]/m_eff$dof),
                         sprintf("$a M_K^\\mathrm{mfit} = %s $",
                                              tex.catwitherror(x=fitresultsrange[x],
                                              dx=fitresultserr[x],
                                              digits=3,
                                              with.dollar=FALSE)),
                         sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                                              fitresultschisq[x]/(t2-t1+1-2)),
                         sprintf("$a M_K^\\mathrm{weighted~median} = %s $",
                                              tex.catwitherror(x=finalres,
                                              dx=sd(fit_bs),
                                              digits=3,
                                              with.dollar=FALSE))
                        )
              )
              pcol <- col2rgb("blue", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              rect(xleft=t1-1,
                   xright=t2-1,
                   ybottom=fitresultsrange[x]-fitresultserr[x],
                   ytop=fitresultsrange[x]+fitresultserr[x],
                   col=pcol,
                   border=NA
                )
              lines(x=c(t1-1,t2-1),
                    y=c(fitresultsrange[x],fitresultsrange[x]),
                    col="blue"
              )
              pcol <- col2rgb("red", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              rect(xleft=tmin-1,
                   xright=tmax-1,
                   ybottom=finalres-sd(fit_bs),
                   ytop=finalres+sd(fit_bs),
                   col=pcol,
                   border=NA
              )
              lines(x=c(t1-1,t2-1),
                    y=c(finalres,finalres),
                    col="red"
              )

              err<-NULL
              for (t in tmin:tmax){ temp<-lapply(X=1:boot.R,FUN=function(rw_idx){newb$cf.tsboot$t[rw_idx,t]/coshfit(t-1,2*L,c(fitresultsbs1[rw_idx,x],fitresultsbs[rw_idx,x]))});temp2 <-do.call(rbind, temp);err<-rbind(err,sd(as.vector(temp2)))}

              if (useCov == TRUE){
                plotwitherror(x = (tmin-1):(tmax-1),
                              y = newb$cf.tsboot$t0[tmin:tmax]/coshfit((tmin-1):(tmax-1),2*L,c(fitresultsrange1[x],fitresultsrange[x])),
                              dy = err,
                              ylab = "C(t)/f(t)",
                              xlab="$t/a$",
                              main=sprintf("$W=%4f~\\delta=%4f~t_1=%d~t_2=%d$", normalizedweight[x],sum(abs(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))-listt))/length(seq(t1-1,t2-1)),t1-1,t2-1)
                    )
              }
              else if (useCov == FALSE){
                plotwitherror(x = (tmin-1):(tmax-1),
                              y = newb$cf.tsboot$t0[tmin:tmax]/coshfit((tmin-1):(tmax-1),2*L,c(fitresultsrange1[x],fitresultsrange[x])),
                              dy = err,
                              ylab = "C(t)/f(t)",
                              xlab="$t/a$",
                              main=sprintf("$Id=%d~\\delta=%4f~t_1=%d~t_2=%d$", which.min(listdev),sum(abs(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))-listt))/length(seq(t1-1,t2-1)),t1-1,t2-1)
                      )

              }
              pcol <- col2rgb("red", alpha = TRUE)/255
              pcol[4] <- 0.35
              pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
              erma<-max(err)
              rect(xleft=t1-1,
                   xright=t2-1,
                   ybottom=min(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))-erma),
                   ytop=max(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))+erma),
                   col=pcol,
                   border=NA
              )

              x <- x + 1
            }#t2
         }#t1
       }#useCov
        tikz.finalize(tikzfiles)
      }#mus
    }#L
  }#mul
}#beta
save(obslistkaon,
     file = "obskaon.RData")
save(kaonmass_summary,
     file = "obskaon_summary.RData")
