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




pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F)

ggplot(pop , aes(x = Surface)) + geom_histogram()

ggplot(pop , aes(sample = Surface)) + geom_qq() + geom_qq_line(col = "red")

load("../donnees/opto")

moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

ggplot(moy_geno , aes(x = Surface)) + geom_histogram()
ggplot(moy_geno , aes(sample = Surface)) + geom_qq() + geom_qq_line(col = "red")

# traits dont on regarde l'évolution
traits <- c("preco","N_flag","hauteur")






# en selectionnant selon un seuil -----------------------------------------


# seuils de troncation
seuils <- seq(19, 22 , 0.5)

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
    selection[j,"nb_geno"] <- length(unique(popsel$geno))
    selection[j,"i"] <- i
    selection[j,"a"] <- a
    j <- j+1
  }
  
}


# Selection sur moyenne d'un genotype

selection_lot <- data.frame()
j <- 1

for (a in seuils){
  
  gensel <- moy_geno %>% filter(Surface > a)
  
  popsel <- pop %>% filter(geno %in% gensel$geno)
  
  p <- nrow(popsel)/nrow(pop)
  i <- i_p(p)
  
  sel <- apply(X = popsel[,traits] , MARGIN = 2 , FUN = calculs , don = popsel , simplify = T)
  
  for (t in traits){
    b <- as.data.frame(sel[t])
    selection_lot[j,"trait"] <- t
    selection_lot[j,"H2"] <- b[1,paste0(t,".H2")]
    selection_lot[j,"R"] <- b[1,paste0(t,".mu")] - mean(pop[,t] , na.rm = T)
    selection_lot[j,"sigma"] <- b[1,paste0(t,".sigma")]
    selection_lot[j,"vg"] <- b[1,paste0(t,".vg")]
    selection_lot[j,"nb_geno"] <- length(unique(popsel$geno))
    selection_lot[j,"i"] <- i
    selection_lot[j,"a"] <- a
    j <- j+1
  }
  
}


# resultats 

selection$type <- "grain"
selection_lot$type <- "lot"

res <- rbind(selection , selection_lot)

ggplot(res , aes(x = nb_geno , fill = type , col = type)) + geom_histogram()

ggplot(res , aes(x = a , y = i , col = type)) + geom_point() + geom_line()
# Intensités de sélection plus forte en lot car on n'a pas planté assez de plantes

res2 <- gather(res  , key = "variable" , value = "value" , 2:5)

ggplot(res2 %>% filter(variable == "H2") , aes(x = i , y = value)) + geom_line() + geom_point() + facet_wrap(~trait)

ggplot(res %>% filter(trait == "hauteur") , aes(x = i , y = vg , col = type)) + geom_point() + geom_line()
ggplot(res %>% filter(trait == "hauteur") , aes(x = i , y = sigma , col = type)) + geom_point() + geom_line()
ggplot(res %>% filter(trait == "hauteur") , aes(x = i , y = R , col = type)) + geom_point() + geom_line()

ggplot(res %>% filter(trait == "preco") , aes(x = i , y = vg , col = type)) + geom_point() + geom_line()
ggplot(res %>% filter(trait == "preco") , aes(x = i , y = sigma , col = type)) + geom_point() + geom_line()
ggplot(res %>% filter(trait == "preco") , aes(x = i , y = R , col = type)) + geom_point() + geom_line()








# en selectionnant un nombre donné de grains ------------------------------

nsel <- 300
