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

pred_germi_nirs <- data.frame()

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
  pred_germi_nirs <- rbind(pred_germi_nirs,cv_err)
}

rm(cv_err,seg,sp,X,Y,i)

save(pred_germi_nirs , file = "pred_germi_nirs")


ggplot(pred_germi_nirs , aes(x = nlv , y = y1)) + geom_point() + labs(title = "Taux d'erreur de la prédiction en fonction du nombre de variables latentes" , x = "Nombre de variable latentes" , y = "Taux d'erreur") + facet_wrap(~traitement)



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






# Prediction variance du grain --------------------------------------------
rm(list=ls())

load("../donnees/opto")
load("../donnees/spectres")


lignes <- merge(opto , spectres$brutes , by = "row.names")$Row.names

# Variances a predire
don <- opto %>% group_by(geno) %>% summarise(var_S = sd(Surface)) %>% column_to_rownames(var = "geno")


# calcul des spectres moyens
spectres_moyens <- list()

i <- 1

for (sp in spectres){
  
  sp <- sp[lignes,]
  sp$geno <- sapply(strsplit(row.names(sp), split = "_") , "[",1)
  
  sp_moy <- data.frame()
  
  for (g in unique(sp$geno)){
    un_sp <- sp %>% filter(geno == g) %>% select(!c("geno")) %>% apply(MARGIN = 2 , FUN = mean)
    sp_moy <- rbind(sp_moy , c(un_sp,g))
    names(sp_moy) <- names(sp)
  }
  
  sp_moyen <- sp_moy %>% column_to_rownames(var = "geno") %>% apply(MARGIN = 2 , FUN = as.numeric) %>% as.data.frame()
  row.names(sp_moyen) <- sp_moy$geno
  
  spectres_moyens[[ names(spectres)[i] ]] <- sp_moyen
  
  print(i)
  
  i <- i+1
}


rm(i,g,un_sp,sp_moy,sp,sp_moyen)

save(spectres_moyens , file = "spectres_moyens")


# modelisation

load("spectres_moyens")

lignes <- merge(don , spectres_moyens$brutes , by = "row.names")$Row.names


model <- function(sp , var){
  
  X <- sp[lignes,]
  
  # definition des groupes de la validation croisee
  seg <- segmkf(n = length(var) , K = 5 , type = "random" , nrep = 2)
  
  # Validation croisée
  cv_rmsecv <- gridcvlv(X = X , 
                        Y = var , 
                        segm = seg , 
                        fun = plskern , 
                        nlv = 1:10 , 
                        score = rmsep , 
                        verb = F)$val_rep
  
  cv_r2 <- gridcvlv(X = X , 
                    Y = var , 
                    segm = seg , 
                    fun = plskern , 
                    nlv = 1:10 , 
                    score = cor2 , 
                    verb = F)$val_rep
  
  cv_rmsecv$traitement <- cv_r2$traitement <- names(spectres_moyens)[i]
  
  cv_rmsecv$score <- "RMSECV" ; cv_r2$score <- "R2"
  
  res <<- rbind(res,cv_rmsecv,cv_r2)
  
  i <<- i+1
  
}


i <- 1
res <- data.frame()

lapply(X = spectres_moyens , FUN = model , var = don[lignes,1])

pred_var_grain_nirs_moyen <- res
save(pred_var_grain_nirs_moyen , file = "pred_var_grain_nirs_moyen")


ggplot(data = res %>% filter(score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + facet_wrap(~traitement) + labs(title = "RMSECV en fonction du nombre de variables latentes et du pré-traitement" , x = "Nombre de variables latentes" , y = "RMSECV")

ggplot(data = res %>% filter(score == "R2") , aes(x = nlv , y = y1)) + geom_point() + facet_wrap(~traitement) + labs(title = "R2 en fonction du nombre de variables latentes et du pré-traitement" , x = "Nombre de variables latentes" , y = "R2")

