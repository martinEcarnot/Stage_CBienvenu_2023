rm(list = ls())

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")
load("../donnees/bac")
load("../donnees/genot_pred")

library(prospectr)
library(signal)
library(rrBLUP)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)


pred <- function(genotype , Nrep , Nout , pheno){

  for (rep in 1:Nrep) {
    print(paste("rep =",rep))
    valid <- sample(length(pheno), Nout)
    train <- setdiff(1:length(pheno) , valid)
    
    phenoTrain <- pheno[train]
    genotype_train <- genotype[train,]
    kin <- matA1[train,train]

    gblup_zu <- mixed.solve(y=phenoTrain, K=kin) # kinship based
    
    gblup_xb <- mixed.solve(y=phenoTrain, Z=genotype_train) # "marker" based
    
    # noms pour remplir le tableau de resultats pour modele zu
    nom <- paste0(noms_trait[j],"_zu_",rep)
    
    # Remplissage du tableau pour modele zu
    res[nom,"trait"] <<- noms_trait[j]
    res[nom,"accuracy"] <<- cor(gblup_zu$u[valid], pheno[valid], use="complete.obs")
    res[nom,"model"] <<- "zu"
    
    # noms pour remplir le tableau de resultats pour modele xb
    nom <- paste0(noms_trait[j],"_xb_",rep)
    
    # Remplissage du tableau pour modele xb
    res[nom,"trait"] <<- noms_trait[j]
    res[nom,"accuracy"] <<- cor(gblup_xb$u[valid], pheno[valid], use="complete.obs")
    res[nom,"model"] <<- "xb"
  }
  
  j <<- j+1
  
}



boss <- function(X,Y,nrep,nout){
  # individus en commun
  lignes <<- row.names(Y)
  
  # Noms des variables
  noms_trait <<- names(Y)
  
  res <<- data.frame()
  j <<- 1
  
  # On laisse que les bons individus
  Y <- Y[lignes,]
  X <- X[lignes,]
  
  # compute kinship matrix
  p <- colMeans(X/2)
  q <- 1-p
  
  genot.scaled <- scale(X, center=2*p, scale=sqrt(4*p*q))
  
  matA1 <<- tcrossprod(genot.scaled) / ncol(genot.scaled)    # tcrossprod(X) is equivalent to X %*% t(X)
  # matA1 is your genomic kinship
  
  apply(Y , MARGIN = 2 , FUN = pred ,
        Nrep = nrep , 
        Nout = nout ,
        genotype = X)
}












calcul_BLUP_bac <- function(i , don){
  mod <- lmer(i ~ (1|geno) + BAC + semis, data = don)
  ranef(mod)$geno
}


# Extraction des BLUPs sur plusieurs variables
don <- bac %>% filter(geno != "INCONNU" & N_flag < 4)

traits <- c("nb_epi","poids_epis")


a <- apply(X = don[,traits] , MARGIN = 2 , FUN = calcul_BLUP_bac , don = don)

BLUP <- as.data.frame(a[1])

for (b in 2:length(a)){
  BLUP <- cbind(BLUP , a[b])
}

names(BLUP) <- names(don[,traits])

rm(a,b)



boss(X = genot_pred , Y = BLUP , nrep = 10 , nout = 40)



ggplot(res , aes(x = trait , y = accuracy^2 , fill = model)) + geom_boxplot() + geom_point(aes(col = model))



resultats_genomique2 <- res

load("resultats_genomique")

resultats_genomique <- rbind(resultats_genomique , resultats_genomique2)

#save(resultats_genomique , file ="resultats_genomique")
