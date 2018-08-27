require(parallel)
require(boot)
require(hadron)
source("framework/utils.R")
source("framework/parameters.R")

obslist <- list()
for (beta in betas){
  betastr <- make.betastr(beta)
  obslist[[betastr]] <- list()
  for (mul in muls[[betastr]] ){
    mulstr <- make.mulstr(mul)
    obslist[[betastr]][[mulstr]] <-list()
    for (L in Ls){
      Lstr <- make.Lstr(L)
      new<-try(readtextcf(file=sprintf("/home/pittler/Bonn_postdoc/pik_scattering/pi_charged_p0_%s%03d_%d.dat",latticelist[[match(beta,betas)]],mul*10000,L), T=2*L, sym=TRUE, path="", skip=1, check.t=0, symmetrise=TRUE))
      if( any(class(new) == "try-error" ) ) {
       next
      }
      tikzfiles <- tikz.init(sprintf("quick_analysis%s%03d_%d",latticelist[[match(beta,betas)]],mul*10000,L), width=4.5, height=4.5)

      str(new)     
 
      tmin <- fitrangepion[[betastr]][[mulstr]][[Lstr]][2]
      tmax <- fitrangepion[[betastr]][[mulstr]][[Lstr]][3]
      lmin <- fitrangepion[[betastr]][[mulstr]][[Lstr]][1]
      lmax <- L

      newb<- bootstrap.cf(cf=new,boot.R=boot.R,boot.l=blocklength[[betastr]][[mulstr]][[Lstr]],seed=551,endcorr=TRUE)
      pvaluefit <- vector()
      fitresultsrange <- vector()
      fitresultsrange1 <- vector()
      fitresultsbs <- vector()
      fitresultsbs1 <- vector()
      fitresultserr <- vector()
      fitresultschisq <- vector ()

      for (t1 in tmin:(tmax-lmin)){
        for (t2 in (t1+lmin):tmax){

          if (useCov == TRUE){
             Cinv<-(solve(cov(newb$cf.tsboot$t[,t1:t2])))
          } else if (useCov == FALSE) {
             cd<-lapply(X=t1:t2,FUN=function(rw_idx){sd(newb$cf.tsboot$t[,rw_idx]);})
             fit_bs <- do.call(rbind, cd)
             Cinv<-diag(as.vector(1/fit_bs^2))
          }
      #Cinv<-diag(1/sd(newb$cf.tsboot$t[,17:31])^2)
          pars.init<- initpion[[betastr]][[mulstr]][[Lstr]]
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
            fitresultsrange <- rbind(fitresultsrange,mfitglobalfinal$par[2])
            fitresultsrange1 <- rbind( fitresultsrange1, mfitglobalfinal$par[1])
            fitresultschisq <- rbind(fitresultschisq, mfitglobalfinal$value)
          }
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
          if (useCov == TRUE){
             fitresultsbs1 <- cbind(fitresultsbs1, fit_bs[,1])
             fitresultsbs <- cbind(fitresultsbs, fit_bs[,2])
             fitresultserr <- rbind(fitresultserr, sd(fit_bs[,2]))  
          }
        }#t2
      }#t1
      minstat <- min(fitresultserr)
      weightlist <- 1:length(fitresultserr)
      weightmeff <- ((1-2*(pvaluefit-0.5)^2)*minstat/fitresultserr)^2
      finalres <- weighted.median(x=fitresultsrange,w=weightmeff)
      cd <- lapply(X=1:boot.R,FUN=function(rw_idx){weighted.median(x=fitresultsbs[rw_idx,],w=weightmeff)})
      fit_bs <- do.call(rbind, cd)
      normalizedweight <- (weightmeff/sum(weightmeff))

      obslist[[betastr]][[mulstr]][[Lstr]] <-list()
      obslist[[betastr]][[mulstr]][[Lstr]][["M_pi_corr"]] <-list(val= finalres, dval=sd(fit_bs), fit_bs=fit_bs)

      x <- 1

      for (t1 in tmin:(tmax-lmin)){
        for (t2 in (t1+lmin):tmax){

          efmassfit <- fit.effectivemass(
                 cf=bootstrap.effectivemass(cf=newb,
                                            type='solve'),
                                            t1=t1,
                                            t2=t2,
                                            useCov=useCov
             )
          m_eff <- efmassfit$effmassfit

          plot(efmassfit,
            ylab="$a M_\\pi^\\mathrm{eff}(t)$",
            xlab="$t/a$",
            main=sprintf("$\\beta %1.2f~a \\mu_l = %.4f~L=%d~Cov=%d~t_1=%d~t_2=%d$", beta,as.numeric(mul),L,useCov,t1,t2),
            ylim=c(m_eff$t0[1]-5*m_eff$se,m_eff$t0[1]+5*m_eff$se))
          legend(x="bottomleft",
            bty='n',
            pch=NA,
            lty=1,
            col=c(rep("black",2),rep("blue",2),rep("red",1)),
            legend=c(sprintf("$a M_\\pi^\\mathrm{eff} = %s $",
                            tex.catwitherror(x=m_eff$t0[1],
                                             dx=m_eff$se,
                                             digits=3,
                                             with.dollar=FALSE)),
                     sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                            m_eff$t0[2]/m_eff$dof),
                     sprintf("$a M_\\pi^\\mathrm{mfit} = %s $",
                            tex.catwitherror(x=fitresultsrange[x],
                                             dx=fitresultserr[x],
                                             digits=3,
                                             with.dollar=FALSE)),
                     sprintf("$\\chi^2/\\mathrm{dof} = %.3f $",
                            fitresultschisq[x]/(t2-t1+1-2)),
                     sprintf("$a M_\\pi^\\mathrm{weighted~median} = %s $",
                            tex.catwitherror(x=finalres,
                                             dx=sd(fit_bs),
                                             digits=3,
                                             with.dollar=FALSE))
                    )
          )
          pcol <- col2rgb("blue", alpha = TRUE)/255
          pcol[4] <- 0.35
          pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
          rect(xleft=t1,
               xright=t2,
               ybottom=fitresultsrange[x]-fitresultserr[x],
               ytop=fitresultsrange[x]+fitresultserr[x],
               col=pcol,
               border=NA
          )
          lines(x=c(t1,t2),
                y=c(fitresultsrange[x],fitresultsrange[x]),
                col="blue"
          )
          pcol <- col2rgb("red", alpha = TRUE)/255
          pcol[4] <- 0.35
          pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
          rect(xleft=tmin,
               xright=tmax,
               ybottom=finalres-sd(fit_bs),
               ytop=finalres+sd(fit_bs),
               col=pcol,
               border=NA
          )
          lines(x=c(t1,t2),
                y=c(finalres,finalres),
                col="red"
          )

          err<-NULL
          for (t in tmin:tmax){ temp<-lapply(X=1:2000,FUN=function(rw_idx){newb$cf.tsboot$t[rw_idx,t]/coshfit(t-1,2*L,c(fitresultsbs1[rw_idx,x],fitresultsbs[rw_idx,x]))});temp2 <-do.call(rbind, temp);err<-rbind(err,sd(as.vector(temp2)))}


          plotwitherror(x = tmin:tmax,
                        y = newb$cf.tsboot$t0[tmin:tmax]/coshfit((tmin-1):(tmax-1),2*L,c(fitresultsrange1[x],fitresultsrange[x])),
                        dy = err,
                        ylab = "C(t)/f(t)",
                        xlab="$t/a$",
                        main=sprintf("$Weight=%4f~t_1=%d~t_2=%d$", normalizedweight[x],t1,t2)
                   )
          pcol <- col2rgb("red", alpha = TRUE)/255
          pcol[4] <- 0.35
          pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
          erma<-max(err)
          rect(xleft=t1,
               xright=t2,
               ybottom=min(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))-erma),
               ytop=max(newb$cf.tsboot$t0[t1:t2]/coshfit((t1-1):(t2-1),2*L,c(fitresultsrange1[x],fitresultsrange[x]))+erma),
               col=pcol,
               border=NA
          )


          x <- x + 1
        }#t2
      }#t1
      tikz.finalize(tikzfiles)
    }#L
  }#mul
}#beta
save(obslist,file = "obspion.RData")

