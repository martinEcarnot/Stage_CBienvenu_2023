rm(list=ls())

library(ggplot2)

# parametres --------------------------------------------------------------

mu <- 5
s_g <- 1
s_inter <- 1
s_intra <- 1.5
s_pos <- 1

surf <- 100  #surface initiale
r <- 0.1  #coefficient de réduction de la surface

NG_E <- 50  #nb grains par epi
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


maxi <- function(a){
  x <- a
  y <- a + 0.01
  
  if (is.nan(tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a)) == T){
    x <- y
    y <- y+0.01
  }
  
  while (tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a) < tronc_epi(x = y , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a)){
    x <- y
    y <- y + 0.01
  }
  
  (x+y)/2
}




# Selection grain a grain -------------------------------------------------
p <- r/(NG_E * NE_P)

a <- qnorm(p = 1-p , mean = mu , sd = s_p)

i <- i_p(p)

H2 <- s_g^2 / s_p^2

S <- i * s_p

R <- H2*S




# Selection epi par epi ---------------------------------------------------
a_epi <- qnorm(p = 1-p , mean = mu , sd = s_epi)

S_epi <- maxi(a = a_epi) - mu

i_epi <- S_epi/s_p

H2_epi <- s_g^2 / (s_g^2 + s_inter^2 + s_pos^2 + s_intra^2/NG_E)

R_epi <- H2_epi*S_epi

NE_O <- p * NE_P * d * surf

NE <- NE_P * d * surf

# resultats ---------------------------------------------------------------

# distributions des grains et des épis
x <- seq(mu - 3.5 * s_p , mu + 3.5 * s_p , 0.1)

don <- data.frame(abs = rep(x,2) , 
                  ord =c(dnorm(x = x , mean = mu , s_p) , dnorm(x = x , mean = mu , sd = s_epi)) ,
                  Distribution = c(rep("Tous les grains" , length(x)) , rep("Moyenne des épis" , length(x))))

ggplot(data = don , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = mu) + labs(x = "Phenotype du grain" , y = "Densité de probabilité") + geom_vline(xintercept = a , col = "blue") + geom_vline(xintercept = a_epi , col = "red")


# Distribution des grains sélectionnés
abs_g <- seq(a , mu + 4 * s_p , 0.1)
abs_e <- seq(S_epi + mu - 2 * s_intra , S_epi + mu + 2 * s_intra , 0.1)

don2 <- data.frame(abs = c(abs_g , abs_e) ,
                   ord = c(dnorm(x = abs_g , mean = mu , sd = s_p) / (1-pnorm(q = a , mean = mu , sd = s_p)) , tronc_epi(x = abs_e , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a_epi)) ,
                   Distribution = c(rep("Sélection sur grain" , length(abs_g)) , rep("Sélection sur épi" , length(abs_e))))



ggplot(data = don2 , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = a , col = "blue") + geom_vline(xintercept = a_epi , col = "red")



data.frame(row.names = c("grains" , "epis") , 
           H2 = c(H2,H2_epi),
           i = c(i,i_epi),
           S = c(S,S_epi),
           R = c(R,R_epi))

paste("Nombre d'épis a sélectionner =" , round(NE_O))

paste("Nombre d'épis dans la parcelle =" , round(NE))

paste("Rapport R_epi/R =" , round(R_epi/R , 2))
