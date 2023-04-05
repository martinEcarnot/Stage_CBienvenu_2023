rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(GGally)
library(lme4)
library(lmerTest)
library(rchemo)

load("../donnees/opto")



# Premiers tests pour voir quelle variable on garde -----------------------

# opto <- opto %>% select("geno" , "grain" , "Longueur" , "Longueur.interieure" , "Largeur" , "Perimetre" , "Perimetre.de.Crofton" , "Perimetre.convexe" , "Dimetre.Eq" , "Surface" , "Surface.convexe" , "Finesse" , "Excentricite" , "Compacite" , "Circularite") %>% mutate_at(.vars = "grain" , .funs = as.factor)


# Matrice de corrélation
# ggcorr(opto)
# opto <- opto %>% select(!c("X"))

# ACP
# res.PCA<-PCA(opto,quali.sup=c(1:2),graph=FALSE)
# plot.PCA(res.PCA,choix='var',title="Graphe des variables de l'ACP" , axes = c(1,3))
# plot.PCA(res.PCA,invisible=c('quali','ind.sup'),habillage='Longueur',title="Graphe des individus de l'ACP",label ='none')

# On peut ne garder que Perimetre, Longueur, Surface, Largeur, Finesse




# ACP ---------------------------------------------------------------------

opto <- opto %>% select("geno" , "grain" , "Longueur" , "Largeur" , "Perimetre" , "Surface" , "Finesse")


ggplot(data = opto , aes(sample = Longueur)) + stat_qq() + stat_qq_line()
ggplot(data = opto , aes(sample = Largeur)) + stat_qq() + stat_qq_line()
ggplot(data = opto , aes(sample = Perimetre)) + stat_qq() + stat_qq_line()
ggplot(data = opto , aes(sample = Surface)) + stat_qq() + stat_qq_line()
ggplot(data = opto , aes(sample = Finesse)) + stat_qq() + stat_qq_line()
ggplot(data = opto , aes(sample = Finesse)) + stat_qq() + stat_qq_line()

ggplot(data = opto , aes(x = Longueur)) + geom_density()
ggplot(data = opto , aes(x = Largeur)) + geom_density()
ggplot(data = opto , aes(x = Perimetre)) + geom_density()
ggplot(data = opto , aes(x = Surface)) + geom_density()
ggplot(data = opto , aes(x = Finesse)) + geom_density()
ggplot(data = opto , aes(x = Finesse)) + geom_density()

res.PCA<-PCA(opto,quali.sup=c(1,2),graph=FALSE)
plot.PCA(res.PCA,choix='var',title="Graphe des variables de l'ACP")
plot.PCA(res.PCA,invisible=c('quali','ind.sup'),habillage='Longueur',title="Graphe des individus de l'ACP",label ='none')

rm(res.PCA)

# Visualisation -----------------------------------------------------------

# Histogrammes
ggpairs(opto[,3:ncol(opto)])




# Heritabilites -----------------------------------------------------------

calcul_H2 <- function(i , don){
  mod <- lmer(i ~ (1|geno) , data = don)
  Vg <- VarCorr(mod)$geno[1]
  Vr <- (sigma(mod))^2
  Vg/(Vg + Vr/12)
}

apply(X = opto[,3:7] , MARGIN = 2 , FUN = calcul_H2 , don = opto)




# Predictions NIRS --------------------------------------------------------

load("../donnees/spectres")


# Récupération des lignes communes entre opto et les spectres
lignes <- merge(opto , spectres$brutes , by = "row.names")$Row.names

opto_pls <- opto[lignes,] %>% select(!c("geno" , "grain"))


# Fonctions pour tout faire tourner

model <- function(sp , var){
  
  X <- sp[lignes,]
  
  # definition des groupes de la validation croisee
  seg <- segmkf(n = length(var) , K = 5 , type = "random" , nrep = 2)
  
  # Validation croisée
  cv_rmsecv <- gridcvlv(X = X , 
                        Y = var , 
                        segm = seg , 
                        fun = plskern , 
                        nlv = 1:10 , 
                        score = rmsep , 
                        verb = F)$val_rep
  
  cv_r2 <- gridcvlv(X = X , 
                    Y = var , 
                    segm = seg , 
                    fun = plskern , 
                    nlv = 1:10 , 
                    score = cor2 , 
                    verb = F)$val_rep
  
  cv_rmsecv$traitement <- cv_r2$traitement <- names(spectres)[i]
  
  cv_rmsecv$score <- "RMSECV" ; cv_r2$score <- "R2"
  
  tmp <<- rbind(tmp,cv_rmsecv,cv_r2)
  
  i <<- i+1
  
  if (i > length(spectres)){i <<- 1}
  
}



boucle <- function(i){
  lapply(X = spectres , FUN = model , var = i)
  tmp$var <<- names(opto_pls)[j]
  res <<- rbind(res,tmp)
  tmp <<- data.frame()
  print(names(opto_pls)[j])
  j <<- j+1
}


# Modelisation
i <- 1
j <- 1
res <- data.frame()
tmp <- data.frame()

apply(X = opto_pls , MARGIN = 2 , FUN = boucle)

save(res , file = "models_opto")



rm(i,j,tmp)


ggplot(res %>% filter(var == "Longueur" & score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 pour la Longueur" , x = "nlv" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(var == "Longueur" & score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV pour la Longueur" , x = "nlv" , y = "R2") + facet_wrap(~traitement)



ggplot(res %>% filter(var == "Largeur" & score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 pour la Largeur" , x = "nlv" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(var == "Largeur" & score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV pour la Largeur" , x = "nlv" , y = "R2") + facet_wrap(~traitement)



ggplot(res %>% filter(var == "Perimetre" & score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 pour la Perimetre" , x = "nlv" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(var == "Perimetre" & score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV pour la Perimetre" , x = "nlv" , y = "R2") + facet_wrap(~traitement)



ggplot(res %>% filter(var == "Surface" & score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 pour la Surface" , x = "nlv" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(var == "Surface" & score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV pour la Surface" , x = "nlv" , y = "R2") + facet_wrap(~traitement)



ggplot(res %>% filter(var == "Finesse" & score == "R2") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "R2 pour la Finesse" , x = "nlv" , y = "R2") + facet_wrap(~traitement)

ggplot(res %>% filter(var == "Finesse" & score == "RMSECV") , aes(x = nlv , y = y1)) + geom_point() + labs(title = "RMSECV pour la Finesse" , x = "nlv" , y = "R2") + facet_wrap(~traitement)