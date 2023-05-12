rm(list = ls())

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")
load("../donnees/spectres")
load("../donnees/opto")


library(prospectr)
library(signal)
library(rrBLUP)
library(tidyverse)

# individus en commun

lignes <- merge(opto , spectres$brutes , by = "row.names")$Row.names %>% sort()




# Validation croisée

un_pretrait <- function(sp , Nrep , Nout , pheno , design = NULL , donnees){
  
  # donnees spectrales
  sp <- sp[lignes,]
  sp <- scale(sp, center=T, scale=T)        # scale absorbance at each wavelength (predictor)
  matH <- tcrossprod(sp)/ncol(sp)     # Compute the hyperspectral similarity matrix
  
  print(paste("pre traitement =",names(spectres[i])))
  
  for (rep in 1:Nrep) {
    print(paste("rep =",rep))
    valid <- sample(length(pheno), Nout)
    phenoTrain <- pheno
    phenoTrain[valid] <- NA
    
    hblup <- mixed.solve(y=phenoTrain, K=matH , Z = design)
    
    
    # noms pour remplir le tableau de resultats
    nom <- paste0(names(donnees[j]),"_",names(spectres[i]),"_",rep)
    
    # Remplissage du tableau
    res[nom,"pretraitement"] <<- names(spectres[i])
    res[nom,"trait"] <<- names(donnees[j])
    res[nom,"rep"] <<- rep
    res[nom,"accuracy"] <<- cor(hblup$u[valid], pheno[valid], use="complete.obs")
  }
  
  i <<- i+1
  
}



un_trait <- function(y , nrep , nout , d = NULL , don){
  
  print(paste("trait =",names(don[j])))
  
  i <<- 1
  
  lapply(spectres , FUN = un_pretrait ,
         Nrep = nrep , 
         Nout = nout ,
         design = d ,
         pheno = y , 
         donnees = don)
  
  j <<- j+1
}






phenotypes <- opto[lignes,] %>% select(Perimetre,Surface)

res <- data.frame()
j <- 1

apply(phenotypes , MARGIN = 2 , FUN = un_trait ,
      nrep = 10 , 
      nout = 100 ,
      don = phenotypes)


pred_pheno <- res

pred_pheno 

save(pred_pheno , file = "pred_pheno")

ggplot(pred_pheno , aes(x = pretraitement , y = accuracy)) + geom_boxplot() + facet_wrap(~trait)



