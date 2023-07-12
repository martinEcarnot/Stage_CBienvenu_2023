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





effets_signif <- function(don , traits , semis){
  
  if (semis == T){base <- "~ geno + semis +"}
  
  if (semis == F){base <- "~ geno +"}

  res <<- data.frame()
  
  for (t in traits){
    
    models <- c()
    
    # effets 1 par 1
    models <- c(models , paste(t,base,"BAC"))
    models <- c(models , paste(t,base,"luz"))
    models <- c(models , paste(t,base,"bordure"))
    models <- c(models , paste(t,base,paste0(t,"_voisin")))
    
    # avec luzerne
    models <- c(models , paste(t,base,"luz/BAC"))
    models <- c(models , paste(t,base,"luz/BAC + BAC/bordure"))
    models <- c(models , paste(t,base,"luz/BAC +",paste0(t,"_voisin")))
    models <- c(models , paste(t,base,"luz/BAC + BAC/bordure +",paste0(t,"_voisin")))
    
    # sans luzerne
    models <- c(models , paste(t,base,"BAC + bordure"))
    models <- c(models , paste(t,base,"BAC + bordure +",paste0(t,"_voisin")))
    models <- c(models , paste(t,base,"BAC/bordure"))
    models <- c(models , paste(t,base,"BAC/bordure +",paste0(t,"_voisin")))
    
    for (m in models){
      mod <- lm(as.formula(m) , data = don)
      
      tab <- data.frame(residuals = mod$residuals , fitted = mod$fitted.values)
      
      a <- ggplot(tab , aes(x = fitted , y = residuals)) + geom_point() + geom_smooth(col = "red" , se = F) + labs(title = "Residuals vs Fitted")
      b <- ggplot(tab , aes(sample = residuals)) + geom_qq() + geom_qq_line(col = "red") + labs(title = "Normal-QQ")
      
      c <- ggarrange(a, b , ncol = 1, nrow = 2)
      
      print(annotate_figure(c , top = text_grob(m , face = "bold")))
      
      test <- anova(mod)
      res[m,"geno"] <<- test["geno" , "Pr(>F)"]
      res[m,"luz"] <<- test["luz" , "Pr(>F)"]
      res[m,"semis"] <<- test["semis" , "Pr(>F)"]
      res[m,"BAC"] <<- test["BAC" , "Pr(>F)"]
      res[m,"luz:BAC"] <<- test["luz:BAC" , "Pr(>F)"]
      res[m,paste0(t,"_voisin")] <<- test[paste0(t,"_voisin") , "Pr(>F)"]
      res[m,"R2"] <<- summary(mod)$adj.r.squared
      res[m,"trait"] <<- t
    }
  }

}





H2 <- function(don , traits , semis){
  
  if (semis == T){base <- "~ (1|geno) + semis +"}
  
  if (semis == F){base <- "~ (1|geno) +"}
  
  res <<- data.frame()
  
  for (t in traits){
    
    models <- c()
    
    # effets 1 par 1
    models <- c(models , paste(t,base,"BAC"))
    models <- c(models , paste(t,base,"luz"))
    models <- c(models , paste(t,base,"bordure"))
    models <- c(models , paste(t,base,paste0(t,"_voisin")))
    
    # avec luzerne
    models <- c(models , paste(t,base,"luz/BAC"))
    models <- c(models , paste(t,base,"luz/BAC + BAC/bordure"))
    models <- c(models , paste(t,base,"luz/BAC +",paste0(t,"_voisin")))
    models <- c(models , paste(t,base,"luz/BAC + BAC/bordure +",paste0(t,"_voisin")))
    
    # sans luzerne
    models <- c(models , paste(t,base,"BAC + bordure"))
    models <- c(models , paste(t,base,"BAC + bordure +",paste0(t,"_voisin")))
    models <- c(models , paste(t,base,"BAC/bordure"))
    models <- c(models , paste(t,base,"BAC/bordure +",paste0(t,"_voisin")))
    
    for (m in models){
      mod <- lmer(as.formula(m) , data = don)
      
      Vg <- VarCorr(mod)$geno[1]
      Vr <- (sigma(mod))^2
      
      
      res[m,"H2"] <<- Vg/(Vg + Vr)
      res[m,"trait"] <<- t
    }
  }
  
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








# BAC : quels effets significatifs ----------------------------------------


load("../donnees/bac")

# On réessaye en enlevant les merdouilles

don <- bac %>% filter(geno != "INCONNU" & appel == "present")

traits <- c("preco","hauteur","N_flag","nb_epi","poids_epis")

effets_signif(don = don , traits = traits , semis = T)


# preco : ok R2 = 0.8
# hauteur : ok R2 = 0.36
# N_flag : ok R2 = 0.23
# nb_epi : pas ouf psq pas vraiment quantitatif R2 = 0.16 - 0.17
# poids_epis : manque homoscedasticite R2 = 0.15

ggplot(don , aes(x = nb_epi , y = poids_epis)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)
# d'où la non homoscedasticite
# donc il faut traiter poids_epis à part en prenant en compte nb_epi qui est très drivé par l'envt

don <- bac %>% filter(geno != "INCONNU" & appel == "present")
mod <- lm(poids_epis ~ geno + luz + semis + nb_epi + bordure , data = don)

summary(mod)$adj.r.squared

tab <- data.frame(residuals = mod$residuals , fitted = mod$fitted.values)

a <- ggplot(tab , aes(x = fitted , y = residuals)) + geom_point() + geom_smooth(col = "red" , se = F) + labs(title = "Residuals vs Fitted")
b <- ggplot(tab , aes(sample = residuals)) + geom_qq() + geom_qq_line(col = "red") + labs(title = "Normal-QQ")
ggarrange(a, b , ncol = 1, nrow = 2)

# Hmmm ça reste pas ouf
# bon on va quand même essayer de calculer des BLUPs comme ça hein on verra bieng


res_2 <- res

R2_2 <- res %>% group_by(trait) %>% summarise(moy = mean(R2) , max = max(R2))




# Même chose en enlevant le semis 2

don <- bac %>% filter(geno != "INCONNU" & appel == "present" & semis == "06/01")

traits <- c("preco","hauteur","N_flag","nb_epi","poids_epis")

effets_signif(don = don , traits = traits , semis = F)

R2_3 <- res %>% group_by(trait) %>% summarise(moy = mean(R2) , max = max(R2))

res_3 <- res


# preco : ok R2 = 0.35
# hauteur : ok R2 = 0.3
# N_flag : ok R2 = 0.25
# nb_epi : pas ouf psq pas vraiment quantitatif R2 = 0.24
# poids_epis : manque homoscedasticite R2 = 0.14


R2_2
R2_3


res_22 <- res_2 %>% select(!c("trait","R2"))
r22 <- c()
for (i in 1:nrow(res_22)){
  for (j in 1:ncol(res_22)){
    if (is.na(res_22[i,j]) == F){
      if (res_22[i,j] >= 0.05){
        print(paste("NON /",rownames(res_22)[i],"/",colnames(res_22)[j]))
        r22 <- c(r22 , paste("NON /",rownames(res_22)[i],"/",colnames(res_22)[j]))}
    }
  }
}


res_33 <- res_3 %>% select(!c("trait","R2"))
r33 <- c()
for (i in 1:nrow(res_33)){
  for (j in 1:ncol(res_33)){
    if (is.na(res_33[i,j]) == F){
      if (res_33[i,j] >= 0.05){
        print(paste("NON /",rownames(res_33)[i],"/",colnames(res_33)[j]))
        r33 <- c(r33 , paste("NON /",rownames(res_33)[i],"/",colnames(res_33)[j]))
        }
    }
  }
}






# pas d'effet de la luzerne sur la hauteur
# pas d'effet du bac sur le poids de l'épi
# pas d'effet de covariable spatiale sur le poids de l'epi

# tout le reste est significatif



# Modeles avec les meilleurs R2 :

for (t in traits){
  tab <- res_2 %>% filter(trait == t)
  print(row.names(tab[which(tab$R2 == max(tab$R2)),]))
}


for (t in traits){
  tab <- res_3 %>% filter(trait == t)
  print(row.names(tab[which(tab$R2 == max(tab$R2)),]))
}


don <- bac %>% filter(geno != "INCONNU" & N_flag < 4 & N_flag_voisin < 4 & semis == "06/01")

ggplot(don , aes(x = luz , y = preco)) + geom_boxplot()
ggplot(don , aes(x = luz , y = hauteur)) + geom_boxplot()
ggplot(don , aes(x = luz , y = N_flag)) + geom_boxplot()
ggplot(don , aes(x = luz , y = nb_epi)) + geom_boxplot()
ggplot(don , aes(x = luz , y = poids_epis)) + geom_boxplot()

ggplot(don , aes(x =  BAC , y = preco , fill = luz)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = hauteur , fill = luz)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = N_flag , fill = luz)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = nb_epi , fill = luz)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = poids_epis , fill = luz)) + geom_boxplot()

ggplot(don , aes(x = bordure , y = preco)) + geom_boxplot()
ggplot(don , aes(x = bordure , y = hauteur)) + geom_boxplot()
ggplot(don , aes(x = bordure , y = N_flag)) + geom_boxplot()
ggplot(don , aes(x = bordure , y = nb_epi)) + geom_boxplot()
ggplot(don , aes(x = bordure , y = poids_epis)) + geom_boxplot()

ggplot(don , aes(x =  BAC , y = preco , fill = bordure)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = hauteur , fill = bordure)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = N_flag , fill = bordure)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = nb_epi , fill = bordure)) + geom_boxplot()
ggplot(don , aes(x =  BAC , y = poids_epis , fill = bordure)) + geom_boxplot()

ggplot(don , aes(x = preco_voisin , y = preco)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)
ggplot(don , aes(x = hauteur_voisin , y = hauteur)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)
ggplot(don , aes(x = N_flag_voisin , y = N_flag)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)
ggplot(don , aes(x = nb_epi_voisin , y = nb_epi)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)
ggplot(don , aes(x = poids_epis_voisin , y = poids_epis)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F)

ggplot(don , aes(x = preco_voisin , y = preco)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + facet_wrap(~BAC)
ggplot(don , aes(x = hauteur_voisin , y = hauteur)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + facet_wrap(~BAC)
ggplot(don , aes(x = N_flag_voisin , y = N_flag)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + facet_wrap(~BAC)
ggplot(don , aes(x = nb_epi_voisin , y = nb_epi)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + facet_wrap(~BAC)
ggplot(don , aes(x = poids_epis_voisin , y = poids_epis)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + facet_wrap(~BAC)







# BAC : H2 avec les deux semis --------------------------------------------

don <- bac %>% filter(geno != "INCONNU" & N_flag < 4)

traits <- c("preco","hauteur","N_flag","nb_epi","poids_epis")

H2(don = don , traits = traits , semis = T)







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
