rm(list = ls())

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")

library(rrBLUP)
library(tidyverse)
library(ggplot2)

load("../donnees/spectres")
load("../donnees/spectres_moyens")
load("../donnees/BLUP")
load("../donnees/opto")

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
  Y <- as.data.frame(Y[lignes,])
  
  apply(Y , MARGIN = 2 , FUN = un_trait ,
        nrep = nrep , 
        nout = nout ,
        spectres = X)
}





# tests pour verifier que c'est miex plusieurs spectres sur un BLUP -------



# test pour être sur que ca marche bien

load("../donnees/spectres_moyens")

nrep <- 10
nout <- 40
boss(X = spectres_moy , Y = BLUP %>% select(prot_semis) , nrep = nrep , nout = nout)

ggplot(res , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()
# ok


# avec plein de spectres visant le même BLUP
nrep <- 10
nout <- 400

for (i in 1:nrow(opto)){
  g <- opto[i,"geno"]
  opto[i,"BLUP_prot_semis"] <- BLUP[g,"prot_semis"]
}

don <- opto %>% select(BLUP_prot_semis)

table(sapply(strsplit(row.names(don) , "_") , "[" , 1))
# certains geno avec 6, d'autres avec 12

boss(X = spectres , Y = don , nrep = nrep , nout = nout)

essai <- res

ggplot(essai , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()


# avec un spectre par BLUP

don$geno <- sapply(strsplit(row.names(don) , "_") , "[" , 1)

don$grain <- sapply(strsplit(row.names(don) , "_") , "[" , 2)


don1 <- don %>% filter(grain == 1) %>% select(BLUP_prot_semis)
don2 <- don %>% filter(grain == 2) %>% select(BLUP_prot_semis)
don3 <- don %>% filter(grain == 3) %>% select(BLUP_prot_semis)

nout <- 20

boss(X = spectres , Y = don1 , nrep = nrep , nout = nout)
essai1 <- res


boss(X = spectres , Y = don2 , nrep = nrep , nout = nout)
essai2 <- res


boss(X = spectres , Y = don3 , nrep = nrep , nout = nout)
essai3 <- res


ggplot(essai , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()
ggplot(essai1 , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()
ggplot(essai2 , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()
ggplot(essai3 , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()


# OK bon ça marche mieux avec plusieurs spectres pour prédire le même BLUP










# moyennage des spectres --------------------------------------------------
# 
# # genotypes à conserver
# gen <- row.names(BLUP)
# gen <- gen[-which(gen == 352)] # cf plus loin pour comprendre pourquoi
# 
# 
# # quels genotype ont 6 spectres et lesquels ont 12 ? 
# sp <- spectres$brutes
# sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
# sp$grain <- sapply(strsplit(row.names(sp) , "_") , "[" , 2)
# sp <- sp %>% filter(geno %in% gen)
# 
# 
# a <- table(sp$geno)
# a[which(a<12)]
# # ah non c'est bon on a bien 12 spectres par genotype sauf le 352, on le retirera à chaque fois
# rm(a,sp)
# 
# 
# # fonction pour moyenner les spectres
# 
# moyennage <- function(spectre , nb){
# 
#   spectre$geno <- sapply(strsplit(row.names(spectre) , "_") , "[" , 1)
#   spectre$grain <- sapply(strsplit(row.names(spectre) , "_") , "[" , 2)
#   
#   # initialisation
#   gr <- sample(x = 1:12 , size = nb)
#   tmp <- spectre %>% filter(geno %in% gen & grain %in% gr) %>% mutate(grain = NULL)
#   tmp <- tmp %>% group_by(geno) %>% summarise_all(.funs = mean) 
#   brutes <- tmp
#   
#   # boucle
#   for (i in 1:11){
#     gr <- sample(x = 1:12 , size = nb)
#     tmp <- spectre %>% filter(geno %in% gen & grain %in% gr) %>% mutate(grain = NULL)
#     tmp <- tmp %>% group_by(geno) %>% summarise_all(.funs = mean) 
#     brutes <- rbind(brutes , tmp)
#   }
#   
#   brutes <- brutes %>% arrange(geno) %>% as.data.frame()
#   brutes$grain <- rep(1:12 , nrow(brutes)/12)
#   row.names(brutes) <- paste0(brutes$geno,"_",brutes$grain)
#   
#   brutes %>% mutate(geno = NULL , grain = NULL)
# }
# 
# 
# 
# 
# # spectres moyens par 2
# spectres_moy2 <- lapply(spectres , FUN = moyennage , nb = 2)
# save(spectres_moy2 , file = "../donnees/spectres_moy2")
# rm(spectres_moy2)
# 
# # spectres moyens par 4
# spectres_moy4 <- lapply(spectres , FUN = moyennage , nb = 4)
# save(spectres_moy4 , file = "../donnees/spectres_moy4")
# rm(spectres_moy4)
# 
# # spectres moyens par 6
# spectres_moy6 <- lapply(spectres , FUN = moyennage , nb = 6)
# save(spectres_moy6 , file = "../donnees/spectres_moy6")
# rm(spectres_moy6)
# 
# # spectres moyens par 8
# spectres_moy8 <- lapply(spectres , FUN = moyennage , nb = 8)
# save(spectres_moy8 , file = "../donnees/spectres_moy8")
# rm(spectres_moy8)
# 
# # spectres moyens par 10
# spectres_moy10 <- lapply(spectres , FUN = moyennage , nb = 10)
# save(spectres_moy10 , file = "../donnees/spectres_moy10")
# rm(spectres_moy10)



# prediction --------------------------------------------------------------

# avec moyenne des 12 spectres
nrep <- 10
nout <- 40
boss(X = spectres_moy , Y = BLUP %>% select(prot_semis) , nrep = nrep , nout = nout)

phenomique_grain_lot <- res
phenomique_grain_lot$nb_sp <- "12"

ggplot(phenomique_grain_lot , aes(x = pretraitement , y = accuracy^2)) + geom_boxplot()
# ok

rm(res,spectres_moy,opto)



# creation d'un jeu de donnees de BLUP avec plein de fois le même BLUP pour différents spectres

geno <- row.names(BLUP)[1]
gr <- 1:12
rn <- paste0(geno,"_",gr)
don <- data.frame(row.names = rn , prot_semis = rep(BLUP$prot_semis[1] , 12))

for (i in 2:nrow(BLUP)){
  geno <- row.names(BLUP)[i]
  gr <- 1:12
  rn <- paste0(geno,"_",gr)
  tmp <- data.frame(row.names = rn , prot_semis = rep(BLUP$prot_semis[i] , 12))
  don <- rbind(don,tmp)
}


rm(tmp,rn,gr,BLUP,i,j)

# avec spectres individuels

nrep <- 10
nout <- 400
boss(X = spectres , Y = don , nrep = nrep , nout = nout)

sp1 <- res
sp1$nb_sp <- "1"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp1)

rm(spectres)

save(phenomique_grain_lot , file = "phenomique_grain_lot")


# avec les spectres moyennés


# moyenne de 2 sp
load("../donnees/spectres_moy2")

boss(X = spectres_moy2 , Y = don , nrep = nrep , nout = nout)

sp2 <- res
sp2$nb_sp <- "2"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp2)

rm(spectres_moy2)

save(phenomique_grain_lot , file = "phenomique_grain_lot")

# moyenne de 4 sp
load("../donnees/spectres_moy4")

boss(X = spectres_moy4 , Y = don , nrep = nrep , nout = nout)

sp4 <- res
sp4$nb_sp <- "4"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp4)

rm(spectres_moy4)

save(phenomique_grain_lot , file = "phenomique_grain_lot")

# moyenne de 6 sp
load("../donnees/spectres_moy6")

boss(X = spectres_moy6 , Y = don , nrep = nrep , nout = nout)

sp6 <- res
sp6$nb_sp <- "6"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp6)

rm(spectres_moy6)

save(phenomique_grain_lot , file = "phenomique_grain_lot")

# moyenne de 8 sp
load("../donnees/spectres_moy8")

boss(X = spectres_moy8 , Y = don , nrep = nrep , nout = nout)

sp8 <- res
sp8$nb_sp <- "8"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp8)

rm(spectres_moy8)


save(phenomique_grain_lot , file = "phenomique_grain_lot")

# moyenne de 10 sp
load("../donnees/spectres_moy10")

boss(X = spectres_moy10 , Y = don , nrep = nrep , nout = nout)

sp10 <- res
sp10$nb_sp <- "10"

phenomique_grain_lot <- rbind(phenomique_grain_lot , sp10)

rm(spectres_moy10)


save(phenomique_grain_lot , file = "phenomique_grain_lot")

rm(phenomique_grain_lot)


load("phenomique_grain_lot")


ggplot(phenomique_grain_lot , aes(x = factor(nb_sp , levels = c("1","2","4","6","8","10","12")) , y = accuracy^2)) + geom_boxplot() + facet_wrap(~pretraitement)
