require(parallel)
require(boot)
require(hadron)
source("framework/utils.R")
source("framework/parameters.R")

load("obspionkaon.RData")

make_plots <- FALSE
pikdifference_summary <- data.frame()
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
        if (is.null(obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$val)== TRUE){
          next
        }
        print(obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$val)
        pikdifference_summary <- rbind( pikdifference_summary,
                                    data.frame(
                                           beta=betastr,
                                           mul=mulstr,
                                           L=Lstr,
                                           mus=musstr,
                                           valE1=(obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov0"]]$val-obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$val)/obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$val,
                                           dvalE1=sd((obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov0"]]$fit_bs-obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$fit_bs)/obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E1"]][["cov1"]]$fit_bs),
                                           valE2=(obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov0"]]$val-obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov1"]]$val)/obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov1"]]$val,
                                           dvalE2=sd((obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov0"]]$fit_bs-obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov1"]]$fit_bs)/obslistpionkaon[[betastr]][[mulstr]][[Lstr]][[musstr]][["E2"]][["cov1"]]$fit_bs)
                                    )
                                  )
     }#mus
   }#Lstr
 }#mulstr
}#beta
save(pikdifference_summary,
     file = "pikdifference_summary.RData")


