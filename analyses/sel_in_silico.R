rm(list=ls())

setwd("~/Stage/Analyses")

load("../donnees/bac")

library(tidyverse)
library(lme4)

# fonctions

calculs <- function(i , don){
  mod <- lmer(i ~ (1|geno) + BAC, data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  
  list(H2 = Vg/(Vg + Vr) , vg = Vg , mu = mean(i , na.rm = T) , sigma = sd(i , na.rm = T))
}

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}




pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F) %>% arrange(Surface)
# On met le tableau selon les surfaces croissantes pour que ce soit plus simple après

ggplot(pop , aes(x = Surface)) + geom_histogram()


# traits dont on regarde l'évolution
traits <- c("preco","N_flag")

# seuils de troncation
seuils <- seq(20, 23 , 0.5)





# Selection sur grain individuel

selection <- data.frame()
j <- 1

for (a in seuils){
  popsel <- pop %>% filter(Surface > a)
  
  p <- nrow(popsel)/nrow(pop)
  i <- i_p(p)
  
  sel <- apply(X = popsel[,traits] , MARGIN = 2 , FUN = calculs , don = popsel , simplify = T)
  
  for (t in traits){
    b <- as.data.frame(sel[t])
    selection[j,"trait"] <- t
    selection[j,"H2"] <- b[1,paste0(t,".H2")]
    selection[j,"R"] <- b[1,paste0(t,".mu")] - mean(pop[,t] , na.rm = T)
    selection[j,"sigma"] <- b[1,paste0(t,".sigma")]
    selection[j,"vg"] <- b[1,paste0(t,".vg")]
    selection[j,"i"] <- i
    j <- j+1
  }
  
}


sel2 <- gather(selection  , key = "variable" , value = "value" , 2:5)

ggplot(sel2 , aes(x = i , y = value)) + geom_line() + geom_point() + facet_grid(variable~trait , scales = "free")

# REGARDER AUSSI LE NOMBRE DE GENOTYPES PRESENTS avant pour voir si on fait vraiment de la selection

# Selection sur moyenne d'un genotype (proche selection sur epi)


