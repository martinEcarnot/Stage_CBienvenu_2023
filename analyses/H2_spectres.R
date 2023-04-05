rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(rchemo)
library(nirsextra)
library(lme4)
library(lmerTest)

load("../donnees/spectres")



calcul_H2 <- function(i , don){
  mod <- lmer(i ~ (1|geno) , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}


res <- data.frame()
i <- 1
for (sp in spectres){
  
  sp$geno <- sapply(strsplit(row.names(sp) , split = "_") , "[",1)
  
  H2 <- apply(X = sp[,-ncol(sp)] , MARGIN = 2 , FUN = calcul_H2 , don = sp)
  tmp <- as.data.frame(H2)
  tmp$lambda <- sapply(strsplit(row.names(tmp) , split = "X") , "[",2)
  tmp$traitement <- names(spectres)[i]
  i <- i+1
  
  res <- rbind(res,tmp)
}

res$lambda <- as.numeric(res$lambda)


sp_moyen <- c()
for (sp in spectres){
  sp_moyen <- c(sp_moyen , apply(sp , MARGIN = 2 , FUN = mean))
}

res$sp_moyen <- sp_moyen

save(res , file = "H2_spectres")

ggplot(res , aes(x = sp_moyen , y = H2)) + geom_point() + labs(title = "Héritabilité des longueur d'onde en fonction de l'absorbance du spectre moyen" , x = "Absorbance moyenne" , y = "H2") + facet_wrap(~traitement , scales = "free")

ggplot(res , aes(x = lambda , y = H2)) + geom_line() + labs(title = "Héritabilité des longueurs d'ondes" , x = "Longueur d'onde" , y = "H2") + facet_wrap(~traitement)


