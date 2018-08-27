make.betastr <- function(beta){
  sprintf("beta%.2f",beta)
}
make.mulstr <- function(mul){
  # the ceiling ensures that the conversion takes place correctly, otherwise
  # as.integer converts 10^4*0.0321 to 320, for example...
  sprintf("mul_%05d", as.integer(ceiling(mul*10^5)))
}
make.musstr <- function(mus){
  # the ceiling ensures that the conversion takes place correctly, otherwise
  # as.integer converts 10^4*0.0321 to 320, for example...
  sprintf("mus_%05d", as.integer(ceiling(mus*10^5)))
}
make.Lstr <- function(L){
  sprintf("L_%d", L)
}

make.key <- function(beta,L,mul,mus,method,cname){
  betastr <- make.betastr(beta)
  Lstr <- make.Lstr(L)
  mulstr <- make.mulstr(mul)
  musstr <- make.musstr(mus)
  sprintf("%s_%s_%s_%s_%s_%s", betastr, Lstr,mulstr, musstr,method, cname)
}

r0_over_b <- function(beta){
  1.0 / exp( -1.6805
             -1.7139*(beta-6)
             +0.8155*(beta-6)^2
             -0.6667*(beta-6)^3 )
}

