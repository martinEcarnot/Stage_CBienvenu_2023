rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(rchemo)
library(nirsextra)

load("../donnees/opto")
load("../donnees/spectres")



# Premier test en prenant tous les spectres -------------------------------


# On utilise la surface comme proxy de la masse

# extraction des variance des surfaces
tmp <- opto %>% group_by(geno) %>% summarise(V = sd(Surface)) %>% column_to_rownames(var = "geno")


don <- data.frame()
for (i in row.names(opto)){
  geno <- opto[i,"geno"]
  don[i,"V"] <- tmp[geno,"V"]
}

rm(i,geno,tmp,opto)

# Prédiction avec les différents spectres :

res <- data.frame()

i <- 1
for (sp in spectres){
  
  # Matrices à utiliser
  data <- merge(don , sp , by = "row.names")
  Y <- data$V
  X <- sp[data$Row.names,]
  rm(data)
  
  # definition des groupes de la validation croisee
  seg <- segmkf(n = length(Y) , K = 5 , type = "random" , nrep = 2)
  
  
  # Validation croisée
  cv_rmsecv <- gridcvlv(X = X , 
                        Y = Y , 
                        segm = seg , 
                        fun = plskern , 
                        nlv = 1:20 , 
                        score = rmsep , 
                        verb = F)$val_rep
  
  cv_r2 <- gridcvlv(X = X , 
                    Y = Y , 
                    segm = seg , 
                    fun = plskern , 
                    nlv = 1:20 , 
                    score = cor2 , 
                    verb = F)$val_rep
  
  
  cv_rmsecv$traitement <- cv_r2$traitement <- names(spectres)[i]
  
  cv_rmsecv$score <- "RMSECV" ; cv_r2$score <- "R2"
  
  
  res <- rbind(res,cv_rmsecv,cv_r2)
  
  print(i)
  i <- i+1
}

rm(cv_r2,cv_rmsecv,seg,sp,X,Y,i)

save(res , file = "pred_var_grain_nirs")

ggplot(res %>% filter(score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 de la prédiction en fonction du pré-traitement et du NLV" , x = "Nombre de variables latentes" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV de la prédiction en fonction du pré-traitement et du NLV" , x = "Nombre de variables latentes" , y = "RMSECV") + facet_wrap(~traitement)








# Deuxieme test en prenant les spectres moyens ----------------------------

load("../donnees/opto")

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

save(res , file = "pred_var_grain_nirs_moyen")


ggplot(data = res %>% filter(score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + facet_wrap(~traitement) + labs(title = "RMSECV en fonction du nombre de variables latentes et du pré-traitement" , x = "Nombre de variables latentes" , y = "RMSECV")

ggplot(data = res %>% filter(score == "R2") , aes(x = nlv , y = y1)) + geom_point() + facet_wrap(~traitement) + labs(title = "R2 en fonction du nombre de variables latentes et du pré-traitement" , x = "Nombre de variables latentes" , y = "R2")



