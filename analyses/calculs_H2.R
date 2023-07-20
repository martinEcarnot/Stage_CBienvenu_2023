rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(GGally)
library(lme4)
library(lmerTest)
library(rchemo)
library(ggpubr)




# H2 des donnees optomachine ----------------------------------------------

calcul_H2_opto <- function(i , don){
  mod <- lmer(i ~ (1|geno), data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}


# données simples
load("../donnees/opto")

traits <- c("Longueur" , "Largeur" , "Perimetre" , "Surface" , "Finesse" , "prot_semis")

apply(X = opto[,traits] , MARGIN = 2 , FUN = calcul_H2_opto , don = opto)



# # variance du grain
# opto$vgrain <- ifelse(opto$grain <= 6 , "1" , "bis")
# vg <- opto %>% group_by(geno,vgrain) %>% summarise(vs = sd(Surface))
# 
# # ne peut pas marcher car 1 seul variance par génotype
# calcul_H2(i = vg$vs , don = vg)
# 
# 
# rm(vg,opto)


rm(opto)


# Calculs H2 bacs avec models selectionnes --------------------------------


calcul_H2 <- function(f,don){
  mod <- lmer(f , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}


load("../donnees/bac")

don <- bac %>% filter(geno != "INCONNU" & appel == "present")

# Hauteur
# model selectionne
f <- hauteur ~ (1|geno) + semis + bordure + hauteur_voisin

calcul_H2(f = f , don = don)


# preco
# model selectionne
f <- preco ~ (1|geno) + semis + BAC + bordure + preco_voisin
calcul_H2(f = f , don = don)


# N_flag
# model selec
f <- N_flag ~ semis + BAC + bordure + N_flag_voisin + nb_voisin # pas d'effet genotype, on essaye d'en mettre un quand même

f <- N_flag ~ semis + BAC + bordure + N_flag_voisin + (1|geno)

calcul_H2(f = f , don = don)


# nb_epi
# mod sel
f <- nb_epi ~ semis + BAC + bordure + nb_epi_voisin # pas geno, on essaye quand même pour voir
f <- nb_epi ~ semis + BAC + bordure + nb_epi_voisin + (1|geno) 

calcul_H2(f = f , don = don)



# poids_epis
# mod sel 
f <- poids_epis ~ semis + BAC + nb_epi # pas geno, on essaye quand même pour voir
f <- poids_epis ~ semis + BAC + nb_epi + (1|geno)

calcul_H2(f = f , don = don)







# extraction des BLUPs pour pred pheno ------------------------------------

# Surface
mod <- lmer(Surface ~ (1|geno), data = opto)
BLUP <- ranef(mod)$geno %>% rename("Surface" = "(Intercept)")

# prot_semis
mod <- lmer(prot_semis ~ (1|geno), data = opto)
B <- ranef(mod)$geno %>% rename("prot_semis" = "(Intercept)")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")


extract_BLUP <- function(f,don,var){
  mod <- lmer(f , data = don)
  a <- ranef(mod)$geno
  colnames(a) <- var
  a
}



# hauteur
f <- hauteur ~ (1|geno) + semis + bordure + hauteur_voisin
B <- extract_BLUP(f = f , don = don , var = "hauteur")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")


# preco
f <- preco ~ (1|geno) + semis + BAC + bordure + preco_voisin
B <- extract_BLUP(f = f , don = don , var = "preco")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")


# N_flag
f <- N_flag ~ semis + BAC + bordure + N_flag_voisin + (1|geno)
B <- extract_BLUP(f = f , don = don , var = "N_flag")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")


# poids_epis
f <- poids_epis ~ semis + BAC + nb_epi + (1|geno)
B <- extract_BLUP(f = f , don = don , var = "poids_epis")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")





BLUP$geno <- row.names(BLUP)

BLUP <- BLUP %>% filter(!geno %in% c("anvergur","ixos","obelix")) %>% select(!c("geno"))


save(BLUP , file = "../donnees/BLUP")





# H2 des axes de l'ACP des spectres ---------------------------------------
rm(list = ls())
load("../donnees/spectres")

calcul_H2_pcasp <- function(i , don){
  mod <- lmer(i ~ (1|geno) , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}

fait <- function(X){
  a <- PCA(X , graph = F)
  b <- as.data.frame(a$ind$coord)
  b$geno <- sapply(strsplit(row.names(b) , split = "_") , "[" , 1)
  apply(X = b[,1:5] , MARGIN = 2 , FUN = calcul_H2_pcasp , don = b)
}

H2 <- data.frame()

# spectres brutes

res <- sapply(spectres , FUN = fait)
res <- as.data.frame(res)

H2 <- gather(res)
H2$dim <- as.factor(rep(1:5,5))

ggplot(H2 , aes(x = dim , y = value)) + geom_col() + facet_wrap(~key) + labs(y = "H2" , x = "Dimension de l'ACP")



















# H2 des spectres ---------------------------------------------------------

# load("../donnees/spectres")
# 
# h2_sp <- data.frame()
# i <- 1
# for (sp in spectres){
#   
#   sp$geno <- sapply(strsplit(row.names(sp) , split = "_") , "[",1)
#   
#   H2 <- apply(X = sp[,-ncol(sp)] , MARGIN = 2 , FUN = calcul_H2 , don = sp)
#   tmp <- as.data.frame(H2)
#   tmp$lambda <- sapply(strsplit(row.names(tmp) , split = "X") , "[",2)
#   tmp$traitement <- names(spectres)[i]
#   i <- i+1
#   
#   h2_sp <- rbind(h2_sp,tmp)
# }
# 
# h2_sp$lambda <- as.numeric(h2_sp$lambda)
# 
# 
# sp_moyen <- c()
# for (sp in spectres){
#   sp_moyen <- c(sp_moyen , apply(sp , MARGIN = 2 , FUN = mean))
# }
# 
# h2_sp$sp_moyen <- sp_moyen
# 
# save(h2_sp , file = "H2_spectres")

load("H2_spectres")

ggplot(h2_sp , aes(x = sp_moyen , y = H2)) + geom_point() + labs(title = "Héritabilité des longueur d'onde en fonction de l'absorbance du spectre moyen" , x = "Absorbance moyenne" , y = "H2") + facet_wrap(~traitement , scales = "free")

ggplot(h2_sp , aes(x = lambda , y = H2)) + geom_line() + labs(title = "Héritabilité des longueurs d'ondes" , x = "Longueur d'onde" , y = "H2") + facet_wrap(~traitement)

