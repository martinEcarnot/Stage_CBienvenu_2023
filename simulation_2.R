rm(list=ls())



# parametres --------------------------------------------------------------

mu <- 5
s_g <- 1
s_inter <- 1
s_intra <- 1
s_pos <- 1
seuil <- 11  # seuil de troncation

S <- 10000  #surface initiale
r <- 0.5  #coefficient de réduction de la surface

NG_E <- 70  #nb grains par epi
NE_P <- 3  #nb epi par plante
d <- 300  # densite de plantation






# calculs utiles ----------------------------------------------------------

s_epi <- sqrt(s_g^2 + s_inter^2 + s_pos^2)
s_p <- sqrt(s_intra^2 + s_g^2 + s_inter^2 + s_pos^2)


# fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

erfc <- function(x) {2 * pnorm(x * sqrt(2), lower = FALSE)}

tronc_epi <- function(x , mu , s_epi , s_intra , seuil){
  -((exp(-((mu-x)^2/(2*(s_epi^2+s_intra^2))))*(s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2+s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2*erf((mu*s_intra^2+s_epi^2*x)/(sqrt(2)*s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2))-((mu*s_intra^2+s_epi^2*x)*erf(sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))+((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)*erf(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2))))))/(sqrt(2*pi)*s_epi*s_intra*(s_epi^2+s_intra^2)*(-2+erfc((-seuil+mu)/(sqrt(2)*s_epi)))))
}


maxi <- function(){
  x <- seuil
  y <- seuil + 0.01
  
  if (is.nan(tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = seuil)) == T){
    x <- y
    y <- y+0.01
  }
  
  while (tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = seuil) < tronc_epi(x = y , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = seuil)){
    x <- y
    y <- y + 0.01
  }
  
  (x+y)/2
}




# Selection grain a grain -------------------------------------------------

p <- r/(NG_E * NE_P)

i <- i_p(p)

H2 <- s_g^2 / s_p^2

S <- i * s_p

R <- H2*S




# Selection epi par epi ---------------------------------------------------

NE_O <- 1

S_epi <- maxi()

i_epi <- S/s_p

H2_epi <- s_g^2 / (s_g^2 + s_inter^2 + s_pos^2 + s_intra^2/NG_E)

R_epi <- H2/S

NE_O <- p * NE_P



# resultats ---------------------------------------------------------------

data.frame(row.names = c("grains" , "epis") , 
           H2 = c(H2,H2_epi),
           i = c(i,i_epi),
           S = c(S,S_epi),
           R = c(R,R_epi))

paste("Nombre d'épis a observer =" , NE_O)
