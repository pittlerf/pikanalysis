require(dplyr)
boot.R <- 2000
boot_sim <- "geom"
useCov <- TRUE


coshfit <-function(t,T,par){
   par[1]*(exp(-t*par[2])+exp(-(T-t)*par[2]))
}
fourpointfit_E2 <- function(t,T,mpion,mkaon,par){

   pollutiont <- exp(-mpion *T)*exp(-(mkaon-mpion)*t)+exp(-mkaon *T)*exp((mkaon-mpion)*t)
   pollutiontp1 <- exp(-mpion *T)*exp(-(mkaon-mpion)*(t+1))+exp(-mkaon *T)*exp((mkaon-mpion)*(t+1))
   par[1]*(exp(-par[2]*t)+exp(-par[2]*(T-t))-pollutiont/pollutiontp1*(exp(-par[2]*(t+1))+exp(-par[2]*(T-(t+1)))))
}
fourpointfit_E1 <- function(t,T,mpion,mkaon,par){
#   par[1] <- 59775.7
#   par[2] <- 0.3804524 
#   par[3] <- -0.01161534
   par[1]*(exp(par[2]*t*-1) + exp(-(T-t)*par[2])-exp(mkaon-mpion)*(exp(-(t+1)*par[2]) + exp(-(T-(t+1))*par[2])))+par[3]*exp((mkaon-mpion)*t)
}
#f(x)=a*(exp(-b*x)+exp(-b*(64-x))-1.145298*(exp(-b*(x+1))+exp(-b*(64-(x+1)))))+c*exp(1.145298*t)
chisqfn_4pt <- function(par,t,T,mpion,mkaon,corrdata,Cinv,fitfn){
  # matrix notation
   res <- (corrdata - fitfn(t,T,mpion,mkaon,par))
   sum( t(res) %*% Cinv %*% res )
}

chisqfn <- function(par,t,T, corrdata,Cinv,fitfn){
  # matrix notation
   res <- (corrdata - fitfn(t,T,par))
   sum( t(res) %*% Cinv %*% res )
}

muls <- list()

muls[[make.betastr(1.90)]] <- c(0.0030, 0.0040, 0.0060, 0.0080, 0.0100)
muls[[make.betastr(1.95)]] <- c(0.0035, 0.0055, 0.0085)
muls[[make.betastr(2.10)]] <- c(0.0030, 0.0045)

betas <- c(1.90,1.95,2.10)
blocklength <- list()
blocklength[[make.betastr(1.90)]] <- list()
blocklength[[make.betastr(1.95)]] <- list()
blocklength[[make.betastr(2.10)]] <- list()

blocklength[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]] <- 3
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]] <- 11
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]] <- 3
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]] <- 2
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]] <- 1
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
blocklength[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]] <- 1

blocklength[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
blocklength[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]] <- 2
blocklength[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
blocklength[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]] <- 1
blocklength[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
blocklength[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]] <- 2


blocklength[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
blocklength[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]] <- 2
blocklength[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
blocklength[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]] <- 3

fitrangekaon <- list()
initkaon <- list()
fitrangekaon[[make.betastr(1.90)]] <- list()
initkaon[[make.betastr(1.90)]] <- list()
fitrangekaon[[make.betastr(1.95)]] <- list()
initkaon[[make.betastr(1.95)]] <- list()
fitrangekaon[[make.betastr(2.10)]] <- list()
initkaon[[make.betastr(2.10)]] <- list()

fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.01850)]]  <- list(444.731,0.2294)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.02250)]]  <- list(420.781,0.249)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.02464)]]  <- list(414.68, 0.26)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]] <- c(7,13,32)

fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(229.683,0.237)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(218.344,0.256)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(214.728,0.266)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]] <- c(5,11,24)

initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.01850)]]  <- list(452.342,0.234)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.02250)]]  <- list(432.802,0.254)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.02464)]]  <- list(422.892,0.264)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]] <- c(7,13,32)

fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(233.692,0.245)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(225.13,0.264)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(220.479,0.274)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]] <- c(5,11,24)

initkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(237.228,0.255)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(229.089,0.273)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(226.5,0.283)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]] <- c(5,11,24)

initkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]]  <- list()
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(239.042,0.264)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(231.688,0.282)
initkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(229.568,0.292)
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
fitrangekaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]] <- c(5,11,24)


initkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]]  <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0160)]]  <- list(390.375,0.206)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0186)]]  <- list(377.269,0.219)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0210)]]  <- list(367.29,0.23)
fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]] <- c(7,13,32)


initkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]]  <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0160)]]  <- list(395.162,0.216)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0186)]]  <- list(383.741,0.228)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0210)]]  <- list(375.208,0.239)
fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]] <- c(7,13,32)


fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
fitrangekaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]] <- c(5,11,24)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]]  <- list()
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0160)]]  <- list(196.454,0.232)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0186)]]  <- list(191.704,0.243)
initkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0210)]]  <- list(187.198,0.253)


fitrangekaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
fitrangekaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]] <- c(10,22,38)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]]  <- list()
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0115)]]  <- list(728.802229,0.150299)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0150)]]  <- list(688.1791079,0.1673344)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0180)]]  <- list(663.7698460,0.1807822)


fitrangekaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
fitrangekaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]] <- c( 7,14,32)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]]  <- list()
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0130)]]  <- list(223.127,0.166)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0150)]]  <- list(218.353,0.176)
initkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0180)]]  <- list(213.3712100,0.1887046)


fitrangepionkaon <- list()
initpionkaon <- list()
fitrangepionkaon[[make.betastr(1.90)]] <- list()
initpionkaon[[make.betastr(1.90)]] <- list()
fitrangepionkaon[[make.betastr(1.95)]] <- list()
initpionkaon[[make.betastr(1.95)]] <- list()
fitrangepionkaon[[make.betastr(2.10)]] <- list()
initpionkaon[[make.betastr(2.10)]] <- list()

fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.01850)]]  <- list(296073.145824268,0.356581278809236, -0.0279129439522094)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.02250)]]  <- list(280911.992592878,0.376410420219657, -0.00921115872951857)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]][[make.musstr(0.02464)]]  <- list(274555.253726533,0.386165683812012,-0.00553893719885271)

fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]] <- c(7,16,28)

fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(69511.5978515027,0.388806535100061,-0.167414067308404)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(67451.8388706706,0.408933400986004,-0.0762733549113995)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(65901.9429062652,0.418633416998375,-0.0518989121328139)
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]] <- c(5,14,22)

initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.01850)]]  <- list(282863.108240769,0.379300290069317,-0.0172140963741271)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.02250)]]  <- list(267961.828276629,0.398336272392773,-0.00605846334511552)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]][[make.musstr(0.02464)]]  <- list(261971.70269088, 0.408118615695931,-0.0034338658039141)
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]] <- c(7,16,28)

fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(64397.1160484382,0.424976742873902,-0.0776869609720234)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(62419.8759578779,0.444202941143327,-0.0385129716995383)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(61284.2933807575,0.453640027086758,-0.0269135437986671)
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]] <- c(5,14,22)

initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(62741.5436542284,0.460895619182365,-0.0351555341455421)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(60541.9003330456,0.478783578849752,-0.0194356505311565)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(60231.6812102298,0.488684111889216,-0.013589879785012)
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]] <- c(5,14,22)

initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]]  <- list()
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.01850)]]  <- list(60095.2533623334,0.492042391762305,-0.0168548956994902)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.02250)]]  <- list(58755.6711553531,0.510090841925263,-0.00999403206327262)
initpionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]][[make.musstr(0.02464)]]  <- list(57927.4624806365,0.51914932046399, -0.00735991816276675)
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
fitrangepionkaon[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]] <- c(5,14,22)


initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]]  <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0160)]]  <- list(203912.852313639,0.334444381366294,-0.0673242705640588)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0186)]]  <- list(198345.089637179,0.34777267043283,-0.0331098070179426)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]][[make.musstr(0.0210)]]  <- list(192823.070046003,0.359334510609848,-0.0177306682337927)
fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]] <- c(7,16,28)


initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]]  <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0160)]]  <- list(184489.513916266,0.374405250860193,-0.0248955554058726)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0186)]]  <- list(179612.639467605,0.387039111917637,-0.0124761293096252)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]][[make.musstr(0.0210)]]  <- list(175320.708284483,0.397885978617381,-0.00742280308815573)
fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]] <- c(7,16,28)


fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
fitrangepionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]] <- c(5,14,22)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]]  <- list()
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0160)]]  <- list(40637.0139747975,0.43257431469522,-0.0477210802223439)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0186)]]  <- list(39633.9774900982,0.444238825498817,-0.0353954097200035)
initpionkaon[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]][[make.musstr(0.0210)]]  <- list(38793.4888318572,0.454287538444116,-0.0258994314461389)


fitrangepionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
fitrangepionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]] <- c(10,27,42)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]]  <- list()
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0115)]]  <- list(702667.238632119,0.250264528135097,-0.0432155606198838)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0150)]]  <- list(665045.467183944,0.267298088649477,-0.0107624765167888)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]][[make.musstr(0.0180)]]  <- list(641752.816772691,0.280696497327184,-0.00346738474744414)


fitrangepionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
fitrangepionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]] <- c( 7,16,28)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]]  <- list()
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0130)]]  <- list(58927.9448053434,0.294725281647923,-0.117163104435799)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0150)]]  <- list(56922.4297963296,0.302956212839376,-0.0804024815217725)
initpionkaon[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]][[make.musstr(0.0180)]]  <- list(54617.3383930194,0.315146124433609,-0.0450602981109212)

initpion <- list()
fitrangepion <- list()
fitrangepion[[make.betastr(1.90)]] <- list()
initpion[[make.betastr(1.90)]] <- list()
fitrangepion[[make.betastr(1.95)]] <- list()
initpion[[make.betastr(1.95)]] <- list()
fitrangepion[[make.betastr(2.10)]] <- list()
initpion[[make.betastr(2.10)]] <- list()


fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
initpion[[make.betastr(1.90)]][[make.mulstr(0.0030)]] <- list()
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0030)]][[make.Lstr(32)]] <- c(7,13,32)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0030)]] [[make.Lstr(32)]] <- c(679.896,0.12288084)

fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
initpion[[make.betastr(1.90)]][[make.mulstr(0.0040)]] <- list()
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]] <- c(5,11,24)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(24)]] <- c(316.434,0.14176206)
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]] <- c(7,13,32)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0040)]][[make.Lstr(32)]] <- c(629,0.14088490)
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
initpion[[make.betastr(1.90)]][[make.mulstr(0.0060)]] <- list()
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]] <- c(5,11,24)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0060)]][[make.Lstr(24)]] <- c(280,0.17128913)


fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
initpion[[make.betastr(1.90)]][[make.mulstr(0.0080)]] <- list()
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]] <- c(5,11,24)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0080)]][[make.Lstr(24)]] <- c(271,0.19820657)
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
initpion[[make.betastr(1.90)]][[make.mulstr(0.0100)]] <- list()
fitrangepion[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]] <- c(5,11,24)
initpion[[make.betastr(1.90)]][[make.mulstr(0.0100)]][[make.Lstr(24)]] <- c(260.03,0.22149677)

fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
initpion[[make.betastr(1.95)]][[make.mulstr(0.0035)]] <- list()
fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]] <- c(7,13,32)
initpion[[make.betastr(1.95)]][[make.mulstr(0.0035)]][[make.Lstr(32)]] <- c(539.83,0.12404299)
fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
initpion[[make.betastr(1.95)]][[make.mulstr(0.0055)]] <- list()
fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]] <- c(7,13,32)
initpion[[make.betastr(1.95)]][[make.mulstr(0.0055)]][[make.Lstr(32)]] <- c(477.74,0.15501990)
fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
initpion[[make.betastr(1.95)]][[make.mulstr(0.0085)]] <- list()
fitrangepion[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]] <- c(5,11,24)
initpion[[make.betastr(1.95)]][[make.mulstr(0.0085)]][[make.Lstr(24)]] <- c(212.355,0.19165412)

fitrangepion[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
initpion[[make.betastr(2.10)]][[make.mulstr(0.0030)]] <- list()
fitrangepion[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]] <- c(10,22,38)
initpion[[make.betastr(2.10)]][[make.mulstr(0.0030)]][[make.Lstr(48)]] <- c(967.3504,0.09737293)
initpion[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
fitrangepion[[make.betastr(2.10)]][[make.mulstr(0.0045)]] <- list()
fitrangepion[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]] <- c( 7,14,32)
initpion[[make.betastr(2.10)]][[make.mulstr(0.0045)]][[make.Lstr(32)]] <- c(272.438,0.12010091)



latticelist <- c("A","B","D")



obsdefs <- list()

obsdefs[["M_K_FSE"]]<- list( name ="M_K_FSE", ltxlabel ="$aM_{K}$") 
obsdefs[["M_pi_FSE"]]<- list( name ="M_pi_FSE", ltxlabel ="$aM_{\\pi}$")
obsdefs[["M_eta"]]<- list( name ="M_eta", ltxlabel ="$aM_{\\eta}$")
obsdefs[["f_pi"]]<- list( name ="f_pi", ltxlabel ="$af_\\text{\\pi}^{\\text{FSE}}$")
obsdefs[["mu_piK_a32"]]<- list( name ="mu_piK_a32", ltxlabel ="$\\mu_{\\pi,K}a_{3/2}$", method=c("E1","E2"))


zpfixinglist <- c ("M1", "M2")
latticespacinglist <- c("a", "r0a")

#r0a[["M1"]] <- list()
#r0a[["M2"]] <- list()

muss <- list()
muss[[make.betastr(1.90)]] <- c(0.01850, 0.02250, 0.02464)
muss[[make.betastr(1.95)]] <- c(0.01600, 0.01860, 0.02100)
muss[[make.betastr(2.10)]] <- c(0.01150, 0.01300, 0.01500, 0.0180)

mussindex <- list()
mussindex[[make.betastr(1.90)]] <- c(1850, 2250, 2464)
mussindex[[make.betastr(1.95)]] <- c(1600, 1860, 2100)
mussindex[[make.betastr(2.10)]] <- c(1150, 1300, 1500, 1800)

methods <- c("E1", "E3" )
Ls <- c(24,32,48)
Ts <- c(48,64,96)
Vs <- outer(X=Ls,Y=Ts,FUN=function(X,Y){ sprintf("L%dT%d",X,Y) } )

# build a data frame of possible parameters
params <- NULL
for( beta in betas ){
  betastr <- make.betastr(beta)
  muss.muls <- expand.grid(muss=muss[[betastr]], muls=muls[[betastr]])
  for( method in methods ){
    for( V in Vs ){
      params <- rbind(params,cbind(muss.muls,data.frame(beta=beta,method=method,V=V)))
    }
  }
}
mypriorsdata=read.csv("/home/pittler/ownCloud/Lattice/Notes/pi-K/Rcode_chiralcontinuum_extrapolation/priors.csv")
#str(mypriorsdata)
zprior<-list()
r0physical<-list()
for( method in zpfixinglist ){
  zprior[[method]]<-list()
  for( beta in betas ){
    betastr<-make.betastr(beta)
    newbeta<-as.double(beta)
    filteredprior<-filter(mypriorsdata,beta==newbeta)
    methodstr<- sprintf("zp%s", method)
    methodstre<- sprintf("zp%se",method) 
    val<-filteredprior[[methodstr]]
    dval<-filteredprior[[methodstre]]
    bs<-rnorm(n=1499,mean=val,sd=dval)
    zprior[[method]][[betastr]]<-list(val=val,
                                      dval=dval,
                                      bs=bs)
    
  }#betas
}#method
mulphysical<-list()
val<-3.70 #3.72
dval<-0.17 # 0.13
bs<-rnorm(n=1499,mean=val,sd=dval)
mulphysical[["M1"]]<-list(val=val,dval=dval,bs=bs)
val<-3.70 #3.63
dval<-0.17 #0.12
bs<-rnorm(n=1499,mean=val,sd=dval)
mulphysical[["M2"]]<-list(val=val,dval=dval,bs=bs)
val<-0.470
dval<-0.012
bs<-rnorm(n=1499,mean=val,sd=dval)
r0physical[["M1"]]<-list(val=val,dval=dval,bs=bs)
val<-0.471
dval<-0.011
bs<-rnorm(n=1499,mean=val,sd=dval)
r0physical[["M2"]]<-list(val=val,dval=dval,bs=bs)
val<-494.2
dval<-0.3
bs<-rnorm(n=1499,mean=val,sd=dval)
kaonmassphysical<-list(val=val,dval=dval,bs=bs)
val<-134.8
dval<-0.3
bs<-rnorm(n=1499,mean=val,sd=dval)
pionmassphysical<-list(val=val,dval=dval,bs=bs)
val<-130.50
dval<-0.13
bs<-rnorm(n=1499,mean=val,sd=dval)
piondecayphysical<-list(val=val,dval=dval,bs=bs)


latticespacing<-list()
for( quantity in latticespacinglist ){
  latticespacing[[quantity]]<-list()
  for( beta in betas ){
    betastr<-make.betastr(beta)
    newbeta<-as.double(beta)
    filteredprior<-filter(mypriorsdata,beta==newbeta)
    quantitye<- sprintf("%se",quantity)
    val<-filteredprior[[quantity]]
    dval<-filteredprior[[quantitye]]
    bs<-rnorm(n=1499,mean=val,sd=dval)
    latticespacing[[quantity]][[betastr]]<-list(val=val,dval=dval,bs=bs)
  }#betas
}#method

interp_cov_types <- c ("diagonal","full")
extrap_cov_types <- c ("diagonal","full")
