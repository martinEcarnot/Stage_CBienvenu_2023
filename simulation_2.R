rm(list=ls())



# parametres --------------------------------------------------------------

mu <- 1
s_g <- 1
s_inter <- 1
s_intra <- 1
s_pos <- 1
seuil <- 1
r <- 0.5


# calculs utiles ----------------------------------------------------------

s_epi <- sqrt(s_g^2 + s_inter^2 + s_pos^2)
s_grain <- sqrt(s_intra^2 + s_g^2 + s_inter^2 + s_pos^2)


# fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

erfc <- function(x) {2 * pnorm(x * sqrt(2), lower = FALSE)}

tronc_epi <- function(x , mu , s_epi , s_intra , seuil){
  -((exp(-((mu-x)^2/(2*(s_epi^2+s_intra^2))))*(s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2+s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2*erf((mu*s_intra^2+s_epi^2*x)/(sqrt(2)*s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2))-((mu*s_intra^2+s_epi^2*x)*erf(sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))+((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)*erf(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2))))))/(sqrt(2*pi)*s_epi*s_intra*(s_epi^2+s_intra^2)*(-2+erfc((-seuil+mu)/(sqrt(2)*s_epi)))))
}


maxi <- function(fun ){
  x
  y
  while (fun(x) < fun(y)){
    x <- y
    y <- y + 0.01
  }
  
  (x+y)/2
}

