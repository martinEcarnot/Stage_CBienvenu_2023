rm(list=ls())

library(ggplot2)

# fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

simulation <- function(mu , vg , vinter , vintra , vpos , r , NG_E , NE_P , NE_O){
  
  # calculs utiles
  s_epi <- sqrt(vg + vinter + vpos)
  s_intra <- sqrt(vintra)
  vp <- vintra + vg + vinter + vpos
  s_p <- sqrt(vp)
  p <- r/(NG_E * NE_P)
  
  # Selection grain a grain 
  ig <- i_p(p)
  H2g <- vg / vp
  Sg <- ig * s_p
  Rg <- H2g*Sg
  
  # Selection epi par epi 
  Se <- ig*s_epi
  ie <- Se/s_p
  H2e <- vg / (vg + vinter + vpos + vintra/NG_E)
  Re <- H2e*Se
  NE_S <- NE_O * p
  

  
  # resultats
  data.frame(mu=mu, vg=vg, vinter=vinter , vintra=vintra , vpos=vpos , r=r , NG_E=NG_E , NE_P=NE_P , NE_O=NE_O , p=p ,
             H2g=H2g , ig=ig , Sg=Sg , Rg=Rg ,
             H2e=H2e , ie=ie , Se=Se , Re=Re)
  

}

# round(
#   data.frame(row.names = c("Sélection sur grain" , "Sélection sur épi"),
#              H2 = c(H2g , H2e),
#              i = c(ig , ie),
#              R = c(Rg , Re),
#              S = c(Sg , Se)),
#   2)



# parametres --------------------------------------------------------------


mu <- 0
vg <- 1
vinter <- 1
vintra <- 0.5
vpos <- 1
r <- 1  #coefficient de réduction de la surface

NG_E <- 70  #nb grains par epi
NE_P <- 3  #nb epi par plante
NE_O <- 5000 #nb epis observés


# analytique et tirage donne la même chose ? ------------------------------

an <- c()
tir <- c()

for (vg in seq(0.8 , 1 , 0.1)){
  for (vinter in seq(0.8 , 1 , 0.1)){
    for (vintra in seq(0.8 , 1 , 0.1)){
      for (vpos in seq(0.8 , 1 , 0.1)){
        for (r in seq(0.8 , 1 , 0.1)){
          for(NG_E in 68 : 70){
            for (NE_P in 1:3){
              for (NE_O in c(4000 , 4500 , 5000)){
                
                tab <- simulation(mu = mu , 
                                  vg = vg , 
                                  vinter = vinter , 
                                  vintra = vintra , 
                                  vpos = vpos , 
                                  r = r , 
                                  NG_E = NG_E , 
                                  NE_P = NE_P , 
                                  NE_O = NE_O)
                
                s_epi <- sqrt(vg + vinter + vpos)
                s_intra <- sqrt(vintra)
                vp <- vintra + vg + vinter + vpos
                s_p <- sqrt(vp)
                p <- r/(NG_E * NE_P)
                
                # Grains
                
                grains <- rnorm(n = 1000000 , mean = mu , sd = s_p)
                ag <- qnorm(p = 1-p , mean = mu , sd = s_p)
                sel <- grains[which(grains>ag)]
                Sg2 <- mean(sel) - mu
                
                # Epis
                
                epis <- rnorm(n = round(1000000/NG_E) , mean = mu , sd = s_epi)
                ae <- qnorm(p = 1-p , mean = mu , sd = s_epi)
                epsel <- epis[which(epis > ae)]
                grsel <- c()
                for (e in epsel){
                  grsel <- c(grsel , rnorm(n = NG_E , mean = e , sd = s_intra))
                }
                
                Se2 <- mean(grsel) - mu
                
              
              an <- c(an , tab[1,"Sg"] , tab[1,"Se"])
              tir <- c(tir , Sg2 , Se2)
              
              }
            }
          }
        }
      }
    }
  }
}


tab <- data.frame(analytique = an , tirage = tir)

ggplot(data = tab , aes(x = analytique , y = tirage)) + geom_point() + geom_abline(slope = 1 , intercept = 0 , col = "red")

cor(tir,an) #0.9968019
