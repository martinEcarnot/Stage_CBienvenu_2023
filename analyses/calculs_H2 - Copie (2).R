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


# H2 des donnees optomachine ----------------------------------------------

calcul_H2_opto <- function(i , don){
  mod <- lmer(i ~ (1|geno), data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}


# données simples
load("../donnees/opto")

traits <- c("Longueur" , "Largeur" , "Perimetre" , "Surface" , "Finesse" , "prot_semis.x")

a <- apply(X = opto[,traits] , MARGIN = 2 , FUN = calcul_H2_opto , don = opto)

H2[traits,"H2"] <- a

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

ggplot(bac , aes(x = hauteur)) + geom_histogram()

don <- bac %>% filter(geno != "INCONNU" & appel == "present" & hauteur > 50 & nb_grain > 10 & semis == "06/01")

don$h <- (don$hauteur - don$hauteur_voisin)/don$hauteur
plot(hauteur~h  , data = don)
hist(don$h)
mean(don$h)

# Hauteur
# model selectionne
f <- poids_epis ~ (1|geno) + bordure

mod <-lmer(f , data = don)

H2["hauteur","H2"] <- calcul_H2(f = f , don = don)


# preco
# model selectionne
f <- preco ~ (1|geno)  + BAC + bordure + nb_voisin + preco_voisin

mod <- lmer(f,data=don)

H2["preco","H2"] <- calcul_H2(f = f , don = don)


# N_flag
# model selec
f <- N_flag ~ BAC + bordure + N_flag_voisin + nb_voisin # pas d'effet genotype, on essaye d'en mettre un quand même

f <- N_flag ~ BAC + bordure + N_flag_voisin + (1|geno)

H2["N_flag","H2"] <- calcul_H2(f = f , don = don)


# nb_epi
# mod sel
f <- nb_epi ~ BAC + bordure + nb_epi_voisin # pas geno, on essaye quand même pour voir
f <- nb_epi ~ BAC + bordure + nb_epi_voisin + (1|geno) 

mod <- lmer(f , data=don)

H2["nb_epi","H2"] <- calcul_H2(f = f , don = don)

ggplot(bac , aes(x = X , y = Y , fill = nb_epi)) + geom_tile() + facet_wrap(~BAC)

# poids_epis
# mod sel 
f <- poids_epis ~ BAC + nb_epi # pas geno, on essaye quand même pour voir
f <- poids_epis ~ BAC + (1|geno)

mod <- lmer(f,data=don)

H2["poids_epis","H2"] <- calcul_H2(f = f , don = don)




# PMG
f <- PMG ~ BAC + bordure + nb_epi + (1|geno)
H2["PMG","H2"] <- calcul_H2(f = f , don = don)



# GSV
f <- GSV ~ bordure + nb_epi + (1|geno)
H2["GSV","H2"] <- calcul_H2(f = f , don = don)


# PMG2
f <- PMG2 ~ BAC + bordure + nb_epi + (1|geno)
H2["PMG2","H2"] <- calcul_H2(f = f , don = don)



# GSV2
f <- GSV2 ~ bordure + nb_epi + (1|geno)
H2["GSV2","H2"] <- calcul_H2(f = f , don = don)


# taille des grains
surface_recolte_moy2 ~ (1|geno)  + BAC + bordure + nb_epi
H2["taille grains","H2"] <- calcul_H2(f = f , don = don)


# prot_recolte
f <- prot_recolte ~ BAC  + N_flag_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin + (1|geno)
H2["prot_recolte","H2"] <- calcul_H2(f = f , don = don)


# nb_grain
f <- nb_grain ~ BAC + (1|geno)
H2["nb_grain","H2"] <- calcul_H2(f = f , don = don)


mod <- lmer(f , data = don)

don$geno2 <- sample(don$geno)

mod2 <- lmer(prot_recolte ~ BAC  + N_flag_voisin + nb_epi_voisin + 
               poids_epis_voisin + nb_voisin + (1|geno2) , data = don)

save(H2 , file="../donnees/H2")






don <- bac %>% filter(geno != "INCONNU")

hist(table(bac$geno,bac$appel)[,1])
table(bac$geno,bac$appel)[,1] + table(bac$geno,bac$appel)[,2]
hist(table(bac$geno,bac$appel)[,2])


hist(table(don$geno,don$appel)[,1])
table(don$geno,don$appel)[,1] + table(don$geno,don$appel)[,2]
hist(table(don$geno,don$appel)[,2])





# H2 bac avec spats -------------------------------------------------------

rm(list=ls())
library(SpATS)
load("../donnees/bac")

don <- bac %>% filter(geno != "INCONNU" & appel == "present") %>% mutate_at(.vars = c("BAC","semis","bordure", "luz") , .funs = as.factor)


# variable x et y adaptees
don$X_spats <- ifelse(don$BAC2 == "x09y04" | don$BAC2 == "x10y04", don$X + 19 , 
                      ifelse(don$BAC2 == "x09y03" | don$BAC2 == "x10y03" , don$X + 39 , don$X))

don$Y_spats <- ifelse(don$BAC2 == "x10y03" | don$BAC2 == "x10y04" | don$BAC2 == "x10y05" , don$Y + 19 , don$Y)

# verif
ggplot(don , aes(x = X_spats , y = Y_spats , fill = BAC2 , col = luz)) + geom_tile()
# c'est good


traits <- c("hauteur","prot_recolte","preco","nb_epi","poids_epis","nb_grain","N_flag","surface_recolte_moy","surface_recolte_moy2","GSV","GSV2","PMG","PMG2")

H2 <- data.frame()

for (t in traits){
  
  # restrictions pour certaines variables
  
  don2 <- don
  
  if (t == "hauteur"){don2 <- don %>% filter(hauteur > 50)}
  
  if (t == "prot_recolte"){don2 <- don %>% filter(prot_recolte < 20 & nb_grain > 10) }
  
  if (t == "nb_grain"){don2 <- don %>% filter(nb_grain > 10)}
  
  model_spat <-SpATS(response = t , 
                      genotype = "geno",
                      spatial = ~ SAP(X_spats, Y_spats, nseg = c(10,10), degree = 3, pord = c(2,2)),
                      genotype.as.random=T, 
                      fixed = ~ semis,
                      data = don2)
  
  plot(model_spat)
  
  H2[t,"H2"] <- getHeritability(model_spat)
  H2[t,"vg"] <- model_spat$var.comp["geno"]
  H2[t,"vx"] <- model_spat$var.comp[1]
  H2[t,"vy"] <- model_spat$var.comp[2]
  
  # 
  # model_spat <-SpATS(response = t , 
  #                    genotype = "geno",
  #                    spatial = ~ SAP(X_spats, Y_spats, nseg = c(10,10), degree = 3, pord = c(2,2)),
  #                    genotype.as.random=T, 
  #                    fixed = ~ semis + BAC,
  #                    data = don2)
  # 
  # plot(model_spat)
  # 
  # 
  # model_spat <-SpATS(response = t , 
  #                    genotype = "geno",
  #                    spatial = ~ SAP(X_spats, Y_spats, nseg = c(10,10), degree = 3, pord = c(2,2)),
  #                    genotype.as.random=T, 
  #                    fixed = ~ semis,
  #                    data = don2)
  # 
  # plot(model_spat)
  # 
  # H2[t,"H2_sans_bac"] <- getHeritability(model_spat)
}






# extraction des BLUPs pour pred pheno lmer ------------------------------------

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
f <- hauteur ~ (1|geno)  + bordure + hauteur_voisin
B <- extract_BLUP(f = f , don = don , var = "hauteur")
BLUP <- merge(BLUP , B , by = "row.names" , all = T) %>% column_to_rownames(var = "Row.names")


# preco
f <- preco ~ (1|geno)  + BAC + bordure + preco_voisin
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










# extraction BLUP spats ---------------------------------------------------

rm(list=ls())
library(SpATS)
load("../donnees/bac")

don <- bac %>% filter(geno != "INCONNU" & appel == "present") %>% mutate_at(.vars = c("BAC","semis","bordure", "luz") , .funs = as.factor)


# variable x et y adaptees
don$X_spats <- ifelse(don$BAC2 == "x09y04" | don$BAC2 == "x10y04", don$X + 19 , 
                      ifelse(don$BAC2 == "x09y03" | don$BAC2 == "x10y03" , don$X + 39 , don$X))

don$Y_spats <- ifelse(don$BAC2 == "x10y03" | don$BAC2 == "x10y04" | don$BAC2 == "x10y05" , don$Y + 15 , don$Y)

# verif
ggplot(don , aes(x = X_spats , y = Y_spats , fill = BAC2 , col = luz)) + geom_tile()
# c'est good


traits <- names(don[c(9,12,14:17,19:23,26:43)])


BLUP_spats <- data.frame(row.names = unique(don$geno))

for (t in traits){
  
  # restrictions pour certaines variables
  
  don2 <- don
  
  if (t == "hauteur"){don2 <- don %>% filter(hauteur > 50)}
  
  if (t == "prot_recolte"){don2 <- don %>% filter(prot_recolte < 20 & nb_grain > 10) }
  
  model_spat <-SpATS(response = t , 
                     genotype = "geno",
                     spatial = ~ SAP(X_spats, Y_spats, nseg = c(10,10), degree = 3, pord = c(2,2)),
                     genotype.as.random=T, 
                     fixed = ~ semis + BAC,
                     data = don2)
  
  #plot(model_spat)
  
  Coeff<-as.data.frame(model_spat$coeff)
  coeff_f<-data.frame(geno=rownames(Coeff)[1:(which(rownames(Coeff)=="Intercept")-1)],
                      BLUP=as.numeric(Coeff[1:(which(rownames(Coeff)=="Intercept")-1),]))
  
  coeff_f <- coeff_f %>% column_to_rownames(var = "geno")
  names(coeff_f) <- t
  
  BLUP_spats <- merge(BLUP_spats , coeff_f , by = "row.names") %>% column_to_rownames(var = "Row.names")
}


save(BLUP_spats , file = "../donnees/BLUP_spats")


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




H2_ACP_sp <- res
save(H2_ACP_sp , file = "../donnees/H2_ACP_sp")











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

