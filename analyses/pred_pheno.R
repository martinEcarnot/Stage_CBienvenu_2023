rm(list = ls())

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")
load("../donnees/spectres")
load("../donnees/spectres_moyens")
load("../donnees/opto")
load("../donnees/bac")

library(prospectr)
library(signal)
library(rrBLUP)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)


# fonctions ---------------------------------------------------------------


pred <- function(sp , Nrep , Nout , pheno){
  
  # donnees spectrales
  sp <- sp[lignes,]
  sp <- scale(sp, center=T, scale=T)        # scale absorbance at each wavelength (predictor)
  
  print(paste("pre traitement =",noms_sp[i]))
  
  for (rep in 1:Nrep) {
    print(paste("rep =",rep))
    valid <- sample(length(pheno), Nout)
    train <- setdiff(1:length(pheno) , valid)
    
    phenoTrain <- pheno[train]
    sp_train <- sp[train,]
    
    
    matH <- tcrossprod(sp_train)/ncol(sp_train)     # Compute the hyperspectral similarity matrix
    
    hblup_zu <- mixed.solve(y=phenoTrain, K=matH) # kinship based
    
    hblup_xb <- mixed.solve(y=phenoTrain, Z=sp_train) # "marker" based
    
    # noms pour remplir le tableau de resultats pour modele zu
    nom <- paste0(noms_trait[j],"_",noms_sp[i],"_zu_",rep)
    
    # Remplissage du tableau pour modele zu
    res[nom,"pretraitement"] <<- noms_sp[i]
    res[nom,"trait"] <<- noms_trait[j]
    res[nom,"accuracy"] <<- cor(hblup_zu$u[valid], pheno[valid], use="complete.obs")
    res[nom,"model"] <<- "zu"
    
    # noms pour remplir le tableau de resultats pour modele xb
    nom <- paste0(noms_trait[j],"_",noms_sp[i],"_xb_",rep)
    
    # Remplissage du tableau pour modele xb
    res[nom,"pretraitement"] <<- noms_sp[i]
    res[nom,"trait"] <<- noms_trait[j]
    res[nom,"accuracy"] <<- cor(hblup_xb$u[valid], pheno[valid], use="complete.obs")
    res[nom,"model"] <<- "xb"
  }
  
  i <<- i+1
  
}



un_trait <- function(y , nrep , nout , spectres){
  
  print(paste("trait =", noms_trait[j]))
  
  i <<- 1
  
  lapply(spectres , FUN = pred ,
         Nrep = nrep , 
         Nout = nout ,
         pheno = y)
  
  j <<- j+1
}


boss <- function(X,Y,nrep,nout){
  # individus en commun
  lignes <<- merge(X[1] , Y , by = "row.names")$Row.names %>% sort()
  
  # Noms des pré traitements
  noms_sp <<- names(X)
  
  # Noms des variables
  noms_trait <<- names(Y)
  
  res <<- data.frame()
  j <<- 1
  
  # On laisse que les bons individus dans Y
  Y <- Y[lignes,]
  
  apply(Y , MARGIN = 2 , FUN = un_trait ,
        nrep = nrep , 
        nout = nout ,
        spectres = X)
}


calcul_BLUP <- function(i , don , x1 = NULL , x2 = NULL , x3 = NULL){
  mod <- lmer(i ~ (1|geno), data = don)
  ranef(mod)$geno
}




# donnees optomachine -------------------------------------------------------

don <- opto


# Extraction des BLUPs sur plusieurs variables
a <- apply(X = opto[,c(3,5,6,10)] , MARGIN = 2 , FUN = calcul_BLUP , don = opto)

BLUP <- as.data.frame(a[1])

for (b in 2:length(a)){
  BLUP <- cbind(BLUP , a[b])
}

names(BLUP) <- names(opto[,c(3,5,6,10)])

rm(a,b)

# verif que les BLUPs suivent bien une loi normale :
ggplot(BLUP , aes(sample=Longueur)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=Largeur)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=Perimetre)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=Surface)) + geom_qq() + geom_qq_line(col = "red")

# verif que BLUP et pheno sont corrélés
phen <- don %>% group_by(geno) %>% summarise_at(.vars = c("Longueur","Largeur","Perimetre","Surface") , .funs = mean) %>% column_to_rownames(var = "geno") %>% merge(BLUP , by = "row.names")

ggplot(phen , aes(x = Surface.x , y = Surface.y)) + geom_point()
ggplot(phen , aes(x = Longueur.x , y = Longueur.y)) + geom_point()
ggplot(phen , aes(x = Largeur.x , y = Largeur.y)) + geom_point()
ggplot(phen , aes(x = Perimetre.x , y = Perimetre.y)) + geom_point()

rm(phen)




# compute entre BLUPs et spectres moyens 


boss(X = spectres_moy , Y = BLUP , nrep=3 , nout=20)


ggplot(res , aes(x = pretraitement , y = accuracy , fill = model)) + geom_boxplot() + facet_wrap(~trait)




# compute en ne moyennant pas les spectres et en mettant son BLUP  a chaque genotype

for (i in 1:nrow(opto)){
  g <- opto[i,"geno"]
  opto[i,"BLUP_longueur"] <- BLUP[g,"Longueur"]
  opto[i,"BLUP_largeur"] <- BLUP[g,"Largeur"]
  opto[i,"BLUP_perimetre"] <- BLUP[g,"Perimetre"]
  opto[i,"BLUP_surface"] <- BLUP[g,"Surface"]
}

boss(X = spectres , Y = opto[,c("BLUP_longueur" , "BLUP_largeur" , "BLUP_perimetre" , "BLUP_surface")] , nrep=3 , nout=400)

pred_pheno_blup <- res

save(pred_pheno_blup , file = "pred_pheno_blup")








# pour donnees bac  -------------------------------------------------------

don <- bac #%>% filter(geno != "INCONNU" & N_flag <= 4)

calcul_BLUP_bac <- function(i , don){
  mod <- lmer(i ~ (1|geno) + BAC + semis, data = don)
  ranef(mod)$geno
}


# Extraction des BLUPs sur plusieurs variables
traits <- c("preco","hauteur")
a <- apply(X = don[,traits] , MARGIN = 2 , FUN = calcul_BLUP , don = don)

BLUP <- as.data.frame(a[1])

for (b in 2:length(a)){
  BLUP <- cbind(BLUP , a[b])
}

names(BLUP) <- names(don[,traits])

rm(a,b)

# verif que les BLUPs suivent bien une loi normale :
ggplot(BLUP , aes(sample=preco)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=hauteur)) + geom_qq() + geom_qq_line(col = "red")

# verif que BLUP et pheno sont corrélés
phen <- don %>% group_by(geno) %>% summarise_at(.vars = c("preco","hauteur") , .funs = mean , na.rm = T) %>% column_to_rownames(var = "geno") %>% merge(BLUP , by = "row.names")

ggplot(phen , aes(x = preco.x , y = preco.y)) + geom_point()
ggplot(phen , aes(x = hauteur.x , y = hauteur.y)) + geom_point()

rm(phen)


# compute entre BLUPs et spectres moyens 


boss(X = spectres_moy , Y = BLUP , nrep=3 , nout=20)


ggplot(res , aes(x = pretraitement , y = accuracy , fill = model)) + geom_boxplot() + facet_wrap(~trait)




# compute en ne moyennant pas les spectres et en mettant son BLUP  a chaque genotype

for (i in 1:nrow( don)){
  g <-  don[i,"geno"]
   don[i,"BLUP_hauteur"] <- BLUP[g,"hauteur"]
   don[i,"BLUP_preco"] <- BLUP[g,"preco"]
}

boss(X = spectres , Y = don[,c("BLUP_hauteur" , "BLUP_preco")] , nrep=3 , nout=400)

ggplot(res , aes(x = pretraitement , y = accuracy , fill = model)) + geom_boxplot() + facet_wrap(~trait)
