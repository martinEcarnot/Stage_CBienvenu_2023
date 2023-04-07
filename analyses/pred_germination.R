setwd("~/Stage/Analyses")
library(ggplot2)
library(tidyverse)
library(rchemo)
library(nirsextra)


rm(list=ls())
load("../donnees/bac")
load("../donnees/spectres")



# Prédiction de la germination par NIRS -----------------------------------


# Prédiction avec les différents spectres :

res <- data.frame()

i <- 1

lignes <- merge(bac , spectres$brutes , by = "row.names")$Row.names

for (sp in spectres){
  
  # Matrices à utiliser
  Y <- bac[lignes,"germi"]
  X <- sp[lignes,]
  
  # definition des groupes de la validation croisee
  seg <- segmkf(n = length(Y) , K = 5 , type = "random" , nrep = 2)
  
  
  # Validation croisée
  cv_err <- gridcvlv(X = X , 
                     Y = Y , 
                     segm = seg , 
                     fun = plsrda , 
                     nlv = 1:20 , 
                     score = err , 
                     verb = F)$val_rep
  
  cv_err$traitement <- names(spectres)[i]
  i <- i+1
  res <- rbind(res,cv_err)
}

rm(cv_err,seg,sp,X,Y,i)

save(res , file = "pred_germi_nirs")


ggplot(res , aes(x = nlv , y = y1)) + geom_point() + labs(title = "Taux d'erreur de la prédiction en fonction du nombre de variables latentes" , x = "Nombre de variable latentes" , y = "Taux d'erreur") + facet_wrap(~traitement)



# Stat Bayes pour voir qualité de la prédiction


rm(list=ls())
load("../donnees/bac")
load("../donnees/spectres")


crossval <- function(sp){
  Ytrain <- bac[train,"germi"]
  Xtrain <- sp[train,]
  
  Ytest <- bac[test,"germi"]
  Xtest <- sp[test,]
  
  mod <- plsrda(X = Xtrain , y = Ytrain , nlv = 1)
  
  pred <- predict.Plsrda(mod , X = Xtest)$pred
  
  tmp <- as.data.frame(cbind(Ytest , pred)) ; tmp$traitement <- names(spectres)[j]
  
  res <<- rbind(res , tmp)
  
  j <<- j+1
}

#https://www.mathemathieu.fr/art/articles-maths/64-thm-bayes
# faut faire la cross val soi même en utilisant plsrda à chaque fois, c'est relou

lignes <- merge(bac , spectres$brutes , by = "row.names")$Row.names

bayes <- data.frame()

for (i in 1:10){
  train <- sample(lignes , size = 1000 , replace = F)
  test <- setdiff(lignes , train)
  
  j <- 1
  res <- data.frame()
  
  lapply(X = spectres , FUN = crossval)$dev2
  
  bayes <- rbind(bayes,res)
  
  print(i)
}

names(bayes) <- c("obs","pred","traitement")

rm(i,j,test,train,res)

bilan <- data.frame(row.names = unique(bayes$traitement) , 
                    vrai_pos = rep(NA,5), 
                    vrai_neg = rep(NA,5), 
                    faux_pos = rep(NA,5), 
                    faux_neg = rep(NA,5))

for (t in unique(bayes$traitement)){
  sub <- subset(bayes , traitement == t)
  
  oui_oui <- nrow(subset(sub , obs == "OUI" & pred == "OUI"))
  oui_non <- nrow(subset(sub , obs == "OUI" & pred == "NON"))
  non_non <- nrow(subset(sub , obs == "NON" & pred == "NON"))
  non_oui <- nrow(subset(sub , obs == "NON" & pred == "OUI"))
  
  vrai_pos <- oui_oui / (oui_oui + non_oui)
  vrai_neg <- non_non / (non_non + oui_non)
  
  faux_pos <- non_oui / (oui_oui + non_oui)
  faux_neg <- oui_non / (non_non + oui_non)
  
  bilan[t,] <- c(vrai_pos,vrai_neg,faux_pos,faux_neg)
}







# Prédiction de la germination par optomachine ----------------------------

rm(list=ls())
load("../donnees/bac")
load("../donnees/opto")

opto <- opto %>% select(!c("geno","grain"))





# Matrices à utiliser
lignes <- merge(bac , opto , by = "row.names")$Row.names
Y <- bac[lignes,"germi"]
X <- opto[lignes,]

# definition des groupes de la validation croisee
seg <- segmkf(n = length(Y) , K = 5 , type = "random" , nrep = 2)


# Validation croisée
res <- gridcvlv(X = X , 
                Y = Y , 
                segm = seg , 
                fun = plsrda , 
                nlv = 1:20 , 
                score = err , 
                verb = F)$val_rep



rm(seg,X,Y)



ggplot(res , aes(x = nlv , y = y1)) + geom_point() + labs(title = "Taux d'erreur de la prédiction en fonction du nombre de variables latentes" , x = "Nombre de variable latentes" , y = "Taux d'erreur")



# Bayes pour voir qualité de la pred

#https://www.mathemathieu.fr/art/articles-maths/64-thm-bayes
# faut faire la cross val soi même en utilisant plsrda à chaque fois, c'est relou

res <- data.frame()

for (i in 1:20){
  train <- sample(lignes , size = 1000 , replace = F)
  test <- setdiff(lignes , train)
  
  Ytrain <- bac[train,"germi"]
  Xtrain <- opto[train,]
  
  Ytest <- bac[test,"germi"]
  Xtest <- opto[test,]
  
  mod <- plsrda(X = Xtrain , y = Ytrain , nlv = 1)
  
  pred <- predict.Plsrda(mod , X = Xtest)$pred
  
  tmp <- as.data.frame(cbind(Ytest , pred))
  
  res <- rbind(res , tmp)
}

names(res) <- c("obs","pred")

unique(res$pred)

rm(i,test,train)

bilan <- data.frame(vrai_pos = NA, 
                    vrai_neg = NA, 
                    faux_pos = NA, 
                    faux_neg = NA)




oui_oui <- nrow(subset(res , obs == "OUI" & pred == "OUI"))
oui_non <- nrow(subset(res , obs == "OUI" & pred == "NON"))
non_non <- nrow(subset(res , obs == "NON" & pred == "NON"))
non_oui <- nrow(subset(res , obs == "NON" & pred == "OUI"))

vrai_pos <- oui_oui / (oui_oui + non_oui)
vrai_neg <- non_non / (non_non + oui_non)

faux_pos <- non_oui / (oui_oui + non_oui)
faux_neg <- oui_non / (non_non + oui_non)

bilan[1,] <- c(vrai_pos,vrai_neg,faux_pos,faux_neg)


