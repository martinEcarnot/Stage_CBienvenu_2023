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
  mod <- lmer(i ~ (1|geno), data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}

# H2 des donnees optomachine ----------------------------------------------

# données simples
load("../donnees/opto")
opto <- opto %>% select("geno" , "grain" , "Longueur" , "Largeur" , "Perimetre" , "Surface" , "Finesse")

apply(X = opto[,3:7] , MARGIN = 2 , FUN = calcul_H2 , don = opto)



# variance du grain
opto$vgrain <- ifelse(opto$grain <= 6 , "1" , "bis")
vg <- opto %>% group_by(geno,vgrain) %>% summarise(vs = sd(Surface))

# ne peut pas marcher car 1 seul variance par génotype
calcul_H2(i = vg$vs , don = vg)


rm(vg,opto)

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








# donnees bac -------------------------------------------------------------

rm(list = ls())

# H2 en prenant en compte la date de semis
calcul_H2_semis <- function(i , don){
  mod <- lmer(i ~ (1|geno) + BAC + semis, data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}



# H2 sans prendre en compte la date de semis (a utiliser avec les donnees du semis 1 seulement)
calcul_H2 <- function(i , don){
  mod <- lmer(i ~ (1|geno) + BAC, data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr)
}


load("../donnees/bac")

bac <- bac %>% filter(geno != "INCONNU")

traits <- c("epiaison","N_flag","preco","hauteur")

apply(X = bac[,traits] , MARGIN = 2 , FUN = calcul_H2_semis , don = bac)

sans_s2 <- bac %>% filter(semis == "06/01")
apply(X = sans_s2[,traits] , MARGIN = 2 , FUN = calcul_H2 , don = sans_s2)


sans_merdouilles <- bac %>% filter(is.na(hauteur) == F & N_flag < 4)
sans_s2_merdouilles <- sans_merdouilles %>% filter(semis == "06/01")
apply(X = sans_merdouilles[,traits] , MARGIN = 2 , FUN = calcul_H2_semis , don = sans_merdouilles)
apply(X = sans_s2_merdouilles[,traits] , MARGIN = 2 , FUN = calcul_H2 , don = sans_s2_merdouilles)




# prot grain semé ---------------------------------------------------------
load("../donnees/opto")

mod <- lmer(prot_semis ~ (1|geno) , data = opto)
Vg <- VarCorr(mod)$geno[1]
Vr <- (sigma(mod))^2
Vg/(Vg + Vr)

BLUP <- ranef(mod)$geno %>% rename(GBLUP = "(Intercept)")

ggplot(BLUP , aes(x = GBLUP)) + geom_histogram()
ggplot(BLUP , aes(sample = GBLUP)) + geom_qq() + geom_qq_line(col = "red")


# Plots -------------------------------------------------------------------


ggplot(bac , aes(x = geno , y = epiaison)) + geom_boxplot()

cherche <- bac %>% group_by(geno) %>% summarise(var_preco = sd(epiaison , na.rm = T))
ggplot(cherche , aes(x = geno , y = var_preco)) + geom_point()

ggplot(bac , aes(x = X , y = Y , fill = epiaison)) + geom_tile() + facet_wrap(~BAC)
ggplot(bac , aes(x = X , y = Y , fill = semis)) + geom_tile() + facet_wrap(~BAC)

for (i in 1:nrow(bac)){
  n <- which(cherche$geno == bac[i,"geno"])
  bac[i,"var_preco"] <- cherche[n,"var_preco"]
}

ggplot(bac , aes(x = X , y = Y , fill = var_preco)) + geom_tile() + facet_wrap(~BAC)


ggplot(bac , aes(x = X , y = Y , fill = appel)) + geom_tile() + facet_wrap(~BAC)

table(bac$appel)/nrow(bac)



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
