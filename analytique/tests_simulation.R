rm(list=ls())

library(ggplot2)

# parametres --------------------------------------------------------------

mu <- 5
vg <- 1.22
vinter <- 0.8
vintra <- 2
vpos <- 1.6
r <- 0.1  #coefficient de réduction de la surface

NG_E <- 50  #nb grains par epi
NE_P <- 3  #nb epi par plante
NE_O <- 5000 #nb epis observés



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
  y <- a + 0.001
  
  if (is.nan(tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a)) == T){
    x <- y
    y <- y+0.001
  }
  
  while (tronc_epi(x = x , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a) < tronc_epi(x = y , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = a)){
    x <- y
    y <- y + 0.01
  }
  
  (x+y)/2
}




# calculs utiles
s_epi <- sqrt(vg + vinter + vpos)
s_intra <- sqrt(vintra)
vp <- vintra + vg + vinter + vpos
s_p <- sqrt(vp)
p <- r/(NG_E * NE_P)

# Selection grain a grain 
ag <- qnorm(p = 1-p , mean = mu , sd = s_p)
ig <- i_p(p)
H2g <- vg / vp
Sg <- ig * s_p
Rg <- H2g*Sg

# Selection epi par epi 
ae <- qnorm(p = 1-p , mean = mu , sd = s_epi)
Se <- maxi(a = ae) - mu
ie <- Se/s_p
H2e <- vg / (vg + vinter + vpos + vintra/NG_E)
Re <- H2e*Se
NE_S <- NE_O * p

res <- round(
  data.frame(row.names = c("Sélection sur grain" , "Sélection sur épi"),
             H2 = c(H2g , H2e),
             i = c(ig , ie),
             R = c(Rg , Re),
             S = c(Sg , Se)),
  2)





# comparaison grosse fonction et juste i ----------------------------------

S <- ig*s_epi
i <- S/s_p
H2 <- H2e
R <- H2 * S

res <- rbind(res , round(c(H2,i,R,S) , 2))

# representations graphiques ----------------------------------------------

# Graph des distributions des grains et des moyennes des épis
x <- seq(mu - 3.5 * s_p , mu + 3.5 * s_p , 0.1)

don <- data.frame(abs = rep(x,2) ,
                  ord =c(dnorm(x = x , mean = mu , sd = s_p) , dnorm(x = x , mean = mu , sd = s_epi)) ,
                  Distribution = c(rep("Grains" , length(x)) , rep("Moyenne des épis" , length(x))))


ggplot(data = don , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = mu) + labs(x = "Phenotype du grain" , y = "Densité de probabilité" , title = "Distributions des grains et des épis") + geom_vline(xintercept = ag , col = "#F8766D") + geom_vline(xintercept = ae , col = "#619CFF") + scale_colour_manual(values = c("#F8766D", "#619CFF")) + theme(plot.background = element_rect(colour = "black"))



## Graph des distributions tronquées

abs_g <- seq(ag , mu + 4 * s_p , 0.1)

abs_e <-   seq(Se + mu - 2 * s_intra , Se + mu + 2 * s_intra , 0.1)

don2 <- data.frame(abs = c(abs_g , abs_e) ,
                   ord = c(dnorm(x = abs_g , mean = mu , sd = s_p) / (1-pnorm(q = ag , mean = mu , sd = s_p)) , tronc_epi(x = abs_e , mu = mu , s_epi = s_epi , s_intra = s_intra , seuil = ae)) ,
                   Distribution = c(rep("Sélection sur grain" , length(abs_g)) , rep("Sélection sur épi" , length(abs_e))))



ggplot(data = don2 , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = ag , col = "#F8766D") + geom_vline(xintercept = ae , col = "#619CFF") + scale_colour_manual(values = c("#619CFF", "#F8766D")) + labs(x = "Phenotype du grain" , y = "Densité de probabilité" , title = "Distributions des grains sélectionnés") + theme(plot.background = element_rect(colour = "black"))




# moyenne des moyennes = moyenne de tout ? --------------------------------


res <- data.frame(moy = 0 , tout = 0)

for(i in 1:20){
  a <- sample(seq(0.1,10,0.1) , 50)
  moy <- c()
  tout <- c()
  for (m in a){
    d <- rnorm(mean = a , sd = s_epi , n=100)
    moy <- c(moy , mean(d))
    tout <- c(tout , d)
  }
res <- rbind(res , c(mean(moy) , mean(tout)))
}

res

unique(res$moy == res$tout)



