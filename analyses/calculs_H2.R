rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(GGally)
library(lme4)
library(lmerTest)
library(rchemo)




# Fonctions ---------------------------------------------------------------

H2_param <- function(i , don , a){
  # fonctionne pour H2 moy
  mod <- lmer(i ~ (1|geno) , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  H2_moy <- Vg/(Vg + Vr/12)
  
  # Knapp et al.
  df_geno <- length(unique(don$geno)) - 1
  df_res <- length(unique(don$geno)) * (12 - 1)
  bas <- 1 - df(x = 1-a/2 , df1 = df_geno , df2 = df_res) * 1/(1-H2_moy)
  haut <- 1 - df(x = a/2 , df1 = df_geno , df2 = df_res) * 1/(1-H2_moy)
  
  print(bas)
  print(H2_moy)
  print(haut)
}

H2_non_param <- function(i , don , a){
  # cette fonction donne l'estimation et l'IC de 1-H2
  opp_H2 <- c()
  
  for (g in unique(don$geno)){
    
    # subset en enlevant un geno à chaque fois
    lignes <- which(don$geno != g)
    j <- i[lignes]
    sub <- don[lignes,] 
    
    # modele
    mod <- lmer(j ~ (1|geno) , data = sub)
    Vg <- VarCorr(mod)$geno[1]
    Vr <- (sigma(mod))^2
    opp_H2 <- c(opp_H2, 1 - Vg/(Vg + Vr))
  }
  
  moy <- mean(opp_H2)
  bas <- moy - dt(x = a , df = length(unique(don$geno)) - 1)
  haut <- moy + dt(x = a , df = length(unique(don$geno)) - 1)
  
  print(bas)
  print(moy)
  print(haut)
  
}

H2_pente <- function(i , don , a){
  mod <- lmer(i ~ (1|geno) , data = don)
  
  g <- ranef(mod)$geno %>% rename("BLUP" = "(Intercept)")
  p <- as.data.frame(cbind(i,don$geno)) %>% mutate_at(.vars = "i" , .funs = as.numeric) %>% group_by(V2) %>% summarise(pheno = mean(i))
  
  mod2 <- lm(g$BLUP ~ p$pheno)
  IC <- confint(mod2 , level = 1-a)
  bas <- IC[2,1]
  haut <- IC[2,2]
  moy <- mod2$coefficients[2]
  
  print(bas)
  print(moy)
  print(haut)
}

calcul_H2 <- function(i , don){
  mod <- lmer(i ~ (1|geno) , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}

# H2 des donnees optomachine ----------------------------------------------

# données simples
load("../donnees/opto")
opto <- opto %>% select("geno" , "grain" , "Longueur" , "Largeur" , "Perimetre" , "Surface" , "Finesse")

apply(X = opto[,3:7] , MARGIN = 2 , FUN = calcul_H2 , don = opto)

rm(opto)



# variance du grain
vgrain <- opto %>% group_by(geno) %>% summarise(vs = sd(Surface))

# ne peut pas marcher car 1 seul variance par génotype
calcul_H2(i = vgrain$vs , don = vgrain)

# H2 des spectres ---------------------------------------------------------

load("../donnees/spectres")

h2_sp <- data.frame()
i <- 1
for (sp in spectres){
  
  sp$geno <- sapply(strsplit(row.names(sp) , split = "_") , "[",1)
  
  H2 <- apply(X = sp[,-ncol(sp)] , MARGIN = 2 , FUN = calcul_H2 , don = sp)
  tmp <- as.data.frame(H2)
  tmp$lambda <- sapply(strsplit(row.names(tmp) , split = "X") , "[",2)
  tmp$traitement <- names(spectres)[i]
  i <- i+1
  
  h2_sp <- rbind(h2_sp,tmp)
}

h2_sp$lambda <- as.numeric(h2_sp$lambda)


sp_moyen <- c()
for (sp in spectres){
  sp_moyen <- c(sp_moyen , apply(sp , MARGIN = 2 , FUN = mean))
}

h2_sp$sp_moyen <- sp_moyen

save(h2_sp , file = "H2_spectres")

ggplot(h2_sp , aes(x = sp_moyen , y = H2)) + geom_point() + labs(title = "Héritabilité des longueur d'onde en fonction de l'absorbance du spectre moyen" , x = "Absorbance moyenne" , y = "H2") + facet_wrap(~traitement , scales = "free")

ggplot(h2_sp , aes(x = lambda , y = H2)) + geom_line() + labs(title = "Héritabilité des longueurs d'ondes" , x = "Longueur d'onde" , y = "H2") + facet_wrap(~traitement)

