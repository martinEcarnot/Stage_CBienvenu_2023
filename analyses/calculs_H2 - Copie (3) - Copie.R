rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(FactoMineR)
#library(GGally)
library(lme4)
library(lmerTest)
#library(rchemo)



H2 <- data.frame()


# Calculs H2 bacs avec models selectionnes --------------------------------


calcul_H2 <- function(f,don){
  mod <- lmer(f , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vb <- VarCorr(mod)$BAC[1]
  Vg/(Vg + Vr + Vb)
}


load("../donnees/bac")


don <- bac %>% filter(geno != "INCONNU" & appel == "present" & hauteur > 50 & nb_grain > 10 & semis == "06/01")


# Hauteur
# model selectionne
f <- poids_epis ~ (1|geno)  + (1|BAC)

mod <-lmer(f , data = don)

H2["hauteur","H2"] <- calcul_H2(f = f , don = don)


# preco
# model selectionne
f <- preco ~ (1|geno)  + (1|BAC)

mod <- lmer(f,data=don)

H2["preco","H2"] <- calcul_H2(f = f , don = don)


# N_flag
# model selec

f <- N_flag ~ (1|geno)  + (1|BAC)

H2["N_flag","H2"] <- calcul_H2(f = f , don = don)


# nb_epi
# mod sel
f <- nb_epi ~ (1|geno)  + (1|BAC)
H2["nb_epi","H2"] <- calcul_H2(f = f , don = don)

# poids_epis
# mod sel 
f <- poids_epis ~ (1|geno)  + (1|BAC)
H2["poids_epis","H2"] <- calcul_H2(f = f , don = don)




# PMG
f <- PMG ~ (1|geno)  + (1|BAC)
H2["PMG","H2"] <- calcul_H2(f = f , don = don)



# GSV
f <- GSV ~ (1|geno)  + (1|BAC)
H2["GSV","H2"] <- calcul_H2(f = f , don = don)


# PMG2
f <- PMG2 ~ (1|geno)  + (1|BAC)
H2["PMG2","H2"] <- calcul_H2(f = f , don = don)



# GSV2
f <- GSV2 ~ (1|geno)  + (1|BAC)
H2["GSV2","H2"] <- calcul_H2(f = f , don = don)


# taille des grains
surface_recolte_moy2 ~ (1|geno)  + (1|BAC)
H2["taille grains","H2"] <- calcul_H2(f = f , don = don)


# prot_recolte
f <- prot_recolte ~ (1|geno)  + (1|BAC)
H2["prot_recolte","H2"] <- calcul_H2(f = f , don = don)


# nb_grain
f <- nb_grain ~ (1|geno)  + (1|BAC)
H2["nb_grain","H2"] <- calcul_H2(f = f , don = don)



