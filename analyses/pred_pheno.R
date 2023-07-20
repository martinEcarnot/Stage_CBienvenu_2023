rm(list = ls())

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")
load("../donnees/spectres")
load("../donnees/spectres_moyens")
load("../donnees/bac")
load("../donnees/BLUP")
load("../donnees/opto")

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
  matH <- tcrossprod(sp)/ncol(sp)     # Compute the hyperspectral similarity matrix
  
  print(paste("pre traitement =",noms_sp[i]))
  
  for (rep in 1:Nrep) {
    print(paste("rep =",rep))
    valid <- sample(length(pheno), Nout)
    
    phenoTrain <- pheno
    phenoTrain[valid] <- NA
    hblup <- mixed.solve(y=phenoTrain, K=matH)
    cor(hblup$u[valid], pheno[valid], use="complete.obs")
    
    # noms pour remplir le tableau de resultats pour modele zu
    nom <- paste0(noms_trait[j],"_",noms_sp[i],"_",rep)
    
    # Remplissage du tableau pour modele zu
    res[nom,"pretraitement"] <<- noms_sp[i]
    res[nom,"trait"] <<- noms_trait[j]
    res[nom,"accuracy"] <<- cor(hblup$u[valid], pheno[valid], use="complete.obs")
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




# pour donnees bac  -------------------------------------------------------

don <- bac %>% filter(geno != "INCONNU" & appel == "present")

traits <- c("hauteur" , "preco" , "N_flag" , "poids_epis")

# verif que les BLUPs suivent bien une loi normale :
ggplot(BLUP , aes(sample=preco)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=hauteur)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=N_flag)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=poids_epis)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=Surface)) + geom_qq() + geom_qq_line(col = "red")
ggplot(BLUP , aes(sample=prot_semis)) + geom_qq() + geom_qq_line(col = "red")

# verif que BLUP et pheno sont corrélés
phen <- don %>% group_by(geno) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T) %>% column_to_rownames(var = "geno") %>% merge(BLUP , by = "row.names")

ggplot(phen , aes(x = preco.x , y = preco.y)) + geom_point()
ggplot(phen , aes(x = hauteur.x , y = hauteur.y)) + geom_point()
ggplot(phen , aes(x = N_flag.x , y = N_flag.y)) + geom_point()
ggplot(phen , aes(x = poids_epis.x , y = poids_epis.y)) + geom_point()

phen <- opto %>% group_by(geno) %>% summarise_at(.vars = c("Surface","prot_semis") , .funs = mean , na.rm = T) %>% column_to_rownames(var = "geno") %>% merge(BLUP , by = "row.names")

ggplot(phen , aes(x = Surface.x , y = Surface.y)) + geom_point()
ggplot(phen , aes(x = prot_semis.x , y = prot_semis.y)) + geom_point()

rm(phen)





# compute entre BLUPs et spectres moyens 

nrep <- 10
nout <- 40
traits <- c("preco" , "hauteur" , "N_flag" , "poids_epis" , "Surface" , "prot_semis")

boss(X = spectres_moy , Y = BLUP , nrep=nrep , nout=nout)


ggplot(res , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot() + facet_wrap(~trait)

blup_sp_moy <- res

blup_sp_moy$donnees <- "blup_sp_moy"








# compute en ne moyennant pas les spectres et en mettant son BLUP  a chaque genotype (ça marche pas du tout)

don <- bac %>% filter(geno != "INCONNU" & appel == "present")

for (i in 1:nrow( don)){
   g <-  don[i,"geno"]
   don[i,"BLUP_N_flag"] <- BLUP[g,"N_flag"]
   don[i,"BLUP_poids_epis"] <- BLUP[g,"poids_epis"]
   don[i,"BLUP_preco"] <- BLUP[g,"preco"]
   don[i,"BLUP_hauteur"] <- BLUP[g,"hauteur"]
}

nout <- 200
traits <- c("BLUP_N_flag","BLUP_poids_epis","BLUP_preco","BLUP_hauteur")

boss(X = spectres , Y = don[,traits] , nrep=nrep , nout=nout)

ggplot(res , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot() + facet_wrap(~trait)

blup_tout_sp <- res
blup_tout_sp$donnees <- "blup_tout_sp"







# Compute sur les donnees pheno directement (marche pas du tout)

nout <- 40
traits <- c("preco" , "hauteur" , "N_flag" , "poids_epis")

boss(X = spectres , Y = don[,traits] , nrep=nrep , nout=nout)

ggplot(res , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot() + facet_wrap(~trait)


phenotype <- res
phenotype$donnees <- "phenotype"




resultats_phenomique <- rbind(blup_sp_moy , blup_tout_sp , phenotype)

#load("resultats_phenomique")

#resultats_phenomique <- rbind(resultats_phenomique , resultats_phenomique2)

save(resultats_phenomique , file = "resultats_phenomique")



























# donnees optomachine -------------------------------------------------------


# compute en ne moyennant pas les spectres et en mettant son BLUP  a chaque genotype

nrep <- 10
nout <- 200

for (i in 1:nrow(opto)){
  g <- opto[i,"geno"]
  opto[i,"BLUP_surface"] <- BLUP[g,"Surface"]
  opto[i,"BLUP_prot_semis"] <- BLUP[g,"prot_semis"]
}

boss(X = spectres , Y = opto[,c("BLUP_surface","BLUP_prot_semis")] , nrep=nrep , nout=nout)

blup_tout_sp <- res
blup_tout_sp$donnees <- "blup_tout_sp"

load("resultats_phenomique")

resultats_phenomique <- rbind(resultats_phenomique,blup_tout_sp)

save(resultats_phenomique , file = "resultats_phenomique")
rm(resultats_phenomique)

# compute avec le phenotype

nout <- 40
traits <- c("prot_semis" , "Surface")

boss(X = spectres , Y = don[,traits] , nrep=nrep , nout=nout)

ggplot(res , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot() + facet_wrap(~trait)


phenotype <- res
phenotype$donnees <- "phenotype"


resultats_phenomique <- rbind(resultats_phenomique, phenotype)

save(resultats_phenomique , file = "resultats_phenomique")
rm(resultats_phenomique)

# resultats ---------------------------------------------------------------

load("resultats_phenomique")


ggplot(resultats_phenomique %>% filter(donnees == "blup_sp_moy") , aes(x = pretraitement , y = accuracy^2)) + facet_wrap(~trait) + geom_boxplot()


ggplot(resultats_phenomique %>% filter(donnees == "blup_tout_sp") , aes(x = pretraitement , y = accuracy^2)) + facet_wrap(~trait) + geom_boxplot()

ggplot(resultats_phenomique %>% filter(donnees == "phenotype") , aes(x = pretraitement , y = accuracy^2)) + facet_wrap(~trait) + geom_boxplot()





