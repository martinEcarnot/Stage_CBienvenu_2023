rm(list=ls())

setwd("~/Stage/Analyses")

load("../donnees/champ")
champ$selection <- relevel(as.factor(champ$selection) , ref = "NT")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(multcomp)
library(lme4)
library(lmerTest)
library(pander)
#library(FactoMineR)



# que est ce qui a ete selectionne ----------------------------------------


load("../donnees/opto_semis_champ")


# qu'est-ce qui a ete selectionne
traits <- c("Surface","Longueur","Largeur","Perimetre","Finesse","Circularite","Excentricite","Index.de.courbure")

NT <- opto_semis_champ %>% filter(taille == "NonTrie")
sel <- opto_semis_champ %>% filter(taille != "NonTrie")

diff <- sel %>% group_by(taille) %>% summarise_at(.vars = traits , .funs = mean) %>% column_to_rownames(var = "taille")

mu <- apply(X = NT[,traits] , MARGIN = 2 , FUN = mean)

sigma <- apply(X = NT[,traits] , MARGIN = 2 , FUN = sd)

i <- data.frame()

for (t in c("Gros","moyen","petit")){
  i[t,traits] <- (diff[t,] - mu)/sigma
}


ggplot(opto_semis_champ , aes(x = Surface)) + geom_histogram()
ggplot(opto_semis_champ , aes(y = Surface , x = taille)) + geom_boxplot()


ggplot(champ , aes(x = surface_recolte_moy)) + geom_histogram()

ggplot(opto_semis_champ , aes(x = Surface , y = Largeur)) + geom_point()
ggplot(opto_semis_champ , aes(x = Surface , y = Longueur)) + geom_point()
ggplot(opto_semis_champ , aes(x = Longueur , y = Largeur)) + geom_point()


ggplot(champ , aes(x = PMG)) + geom_histogram()

ggplot(champ , aes(x = GSV)) + geom_histogram()



# ACP des grains tamises --------------------------------------------------

rm(list=ls())

load("../donnees/opto_semis_champ")

# ça donne rien

# ACP des epis recoltes
load("../donnees/champ")




# gradient spatial ? ------------------------------------------------------

don <- champ %>% filter(!passage %in% c("inconnu1","inconnu2") & !planche %in% c("inconnu1","inconnu2"))


moy <- don %>% group_by(passage,planche) %>% summarise_at(.vars = c("hauteur" , "nb_epillets") , .funs = mean , na.rm = T)


ggplot(moy , aes(x = passage , y = planche , fill = hauteur)) + geom_tile()

ggplot(moy , aes(x = passage , y = planche , fill = nb_epillets)) + geom_tile()


ggplot(champ , aes(x = hauteur , y = nb_epillets)) + geom_point() + geom_smooth(method = "lm" , col = "red")




# effet de la selection ? -------------------------------------------------

rm(list = ls())

load("../donnees/champ")
champ$selection <- relevel(as.factor(champ$selection) , ref = "NT")

ggplot(champ , aes(x = surface_recolte_moy , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = surface_recolte_min , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = surface_recolte_max , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = PMG , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = PMG2 , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = prot_recolte , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = GSV , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = GSV2 , col = selection , fill = selection)) + geom_density(alpha = 0.3)

ggplot(champ , aes(x = nb_grain , y = PMG)) + geom_point() + geom_smooth(method = "lm")
cor(champ$PMG , champ$nb_grain)

ggplot(champ , aes(x = nb_grain)) + geom_histogram()



# avec effet fixe

res <- data.frame()

for (v in names(champ)[c(1,2,7:24)]){
  f <- as.formula(paste(v , "~ nb_grain + selection  + passage"))
  
  mod <- lm(f , data = champ , weights = nb_grain)
  
  # tab <- data.frame(resid=mod$residuals , fit=mod$fitted.values)
  # 
  # print(ggplot(tab , aes(x = fit , y = resid)) + geom_point() + labs(x = "Fitted values" , y = "Residuals" , title = paste(v,"Residuals vs Fitted")))
  # 
  # print(ggplot(tab , aes(sample = resid)) + geom_qq() + geom_qq_line(col = "red") + labs(x = "Theoretical Quantiles" , y = "Residuals" , title = paste(v,"Normal Q-Q")))
  
  
  a <- drop1(mod , .~. , test = "F")["selection" , "Pr(>F)"]
  
  res[v,"pval"] <- a
  res[v,"petit"] <- mod$coefficients["selectionpetit"]
  res[v,"moyen"] <- mod$coefficients["selectionmoyen"]
  res[v,"gros"] <- mod$coefficients["selectiongros"]
}



res %>% filter(pval < 0.05)


# tests de tuckey

library(multcomp)

# PMG
mod <- lm(PMG ~ selection , data = champ)

lettres <- as.data.frame(cld(glht(model = mod , mcp(selection = "Tukey")) , level = 0.05)$mcletters$Letters)

names(lettres)[1] <- "groups"

lettres$selection <- row.names(lettres)

ggplot(champ , aes(x=selection , y=PMG , col=selection, fill=selection))+
  geom_boxplot(alpha = 0.5)+
  theme(legend.position="none") +
  geom_text(data = lettres, aes(label = groups, y = 70 ), size=5) +
  labs(title = "PMG en fonction de la selection" , x = "selection" , y = "PMG") + stat_mean(col = "black")



# PMG2
mod <- lm(PMG2 ~ selection , data = champ)

lettres <- as.data.frame(cld(glht(model = mod , mcp(selection = "Tukey")) , level = 0.05)$mcletters$Letters)

names(lettres)[1] <- "groups"

lettres$selection <- row.names(lettres)

ggplot(champ , aes(x=selection , y=PMG2 , col=selection, fill=selection))+
  geom_boxplot(alpha = 0.5)+
  theme(legend.position="none") +
  geom_text(data = lettres, aes(label = groups, y = 70 ), size=5) +
  labs(title = "PMG2 en fonction de la selection" , x = "selection" , y = "PMG2")




# prot_recolte
mod <- lm(prot_recolte ~ selection , data = champ)

lettres <- as.data.frame(cld(glht(model = mod , mcp(selection = "Tukey")) , level = 0.05)$mcletters$Letters)

names(lettres)[1] <- "groups"

lettres$selection <- row.names(lettres)

ggplot(champ , aes(x=selection , y=prot_recolte , col=selection, fill=selection))+
  geom_boxplot(alpha = 0.5)+
  theme(legend.position="none") +
  geom_text(data = lettres, aes(label = groups, y = 20 ), size=5) +
  labs(title = "prot_recolte en fonction de la selection" , x = "selection" , y = "prot_recolte")

plot(prot_recolte~PMG , data = champ)
plot(prot_recolte~nb_grain , data = champ)
boxplot(GSV~selection , data = champ)




# prot_recolte
mod <- lm(surface_recolte_moy ~ selection , data = champ)

lettres <- as.data.frame(cld(glht(model = mod , mcp(selection = "Tukey")) , level = 0.05)$mcletters$Letters)

names(lettres)[1] <- "groups"

lettres$selection <- row.names(lettres)

ggplot(champ , aes(x=selection , y=surface_recolte_moy , col=selection, fill=selection))+
  geom_boxplot(alpha = 0.5)+
  theme(legend.position="none") +
  geom_text(data = lettres, aes(label = groups, y = 20 ), size=5) +
  labs(title = "prot_recolte en fonction de la selection" , x = "selection" , y = "prot_recolte") + stat_mean(col="black")




# effet sur PMG et prot ! (les anova sont valides, vérifié par les plots)




# on va esayer avec des tests de kruskal pour le reste


res <- data.frame()

for (v in names(champ)[c(1,2,7:24)]){
  f <- as.formula(paste(v , "~ selection"))
  
  mod <- kruskal.test(f , data = champ)
  
  res[v,"pval"] <- mod$p.value
}



res %>% filter(pval < 0.05)



# test de Dunnet pour comparer les selections au temoin

# réalisation de toutes les comparaisons deux à deux sans ajuster les p-values
pp_all <- pairwise.wilcox.test(champ$PMG , champ$selection , p.adjust.method ="none" )

pp_all

# récupération des pvalues des tests de wilcoxon entre la moyenne du groupe contrôle et toutes les autres moyennes.
pp_dunn <- c(pp_all[[3]][,"NT"])
pp_dunn


# ajustement de ces p-values par la méthode de Holm
pp_dunn_adjust <- p.adjust(pp_dunn,method="holm")
pp_dunn_adjust






# avec effet aleatoire

res <- data.frame()

for (v in names(champ)[c(1,2,7:24)]){
  f <- as.formula(paste(v , "~ selection + (1|selection:parcelle)"))
  
  mod <- lmer(f , data = champ)
  
  a <- anova(mod)
  b <- summary(mod)$coefficients
  
  res[v,"pval"] <- a["selection" , "Pr(>F)"]
  res[v,"petit"] <- b["selectionpetit","Estimate"]
  res[v,"moyen"] <- b["selectionmoyen","Estimate"]
  res[v,"gros"] <- b["selectiongros","Estimate"]
}



res %>% filter(pval < 0.05)






# H2 realise --------------------------------------------------------------


champ %>% group_by(parcelle , selection) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T)) %>% filter(selection == "petit")

champ %>% group_by(parcelle , selection) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T)) %>% filter(selection == "moyen")

champ %>% group_by(parcelle , selection) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T)) %>% filter(selection == "gros")

champ %>% group_by(parcelle , selection) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T)) %>% filter(selection == "NT")

# c'est assez homogène entre les parcelles


load("../donnees/S")


# En moyennant tout
moy <- champ %>% group_by(selection) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T) , PMG = mean(PMG , na.rm = T) , PMG2 = mean(PMG2 , na.rm = T)) %>% column_to_rownames(var = "selection")


R <- data.frame()

for (t in c("surface","PMG","PMG2")){
  R["gros",t] <- moy["gros",t] - S["NT",t] 
  R["petit",t] <- moy["petit",t] - S["NT",t]
  R["moyen",t] <- moy["moyen",t] - S["NT",t]
}

S <- S[-4,]

S <- S[sort(row.names(S)),]
R <- R[sort(row.names(R)),]

H2 <- R/S




# parcelle par parcelle

load("../donnees/S")

moy <- champ %>% group_by(selection,parcelle) %>% summarise(surface = mean(surface_recolte_moy , na.rm = T) , PMG = mean(PMG , na.rm = T) , PMG2 = mean(PMG2 , na.rm = T)) %>% mutate(rn = paste0(selection,parcelle)) %>% column_to_rownames(var = "rn")


R2 <- data.frame()
S2 <- data.frame()

for (r in row.names(moy)){
  R2[r,c("surface","PMG","PMG2")] <- moy[r,c("surface","PMG","PMG2")] - S["NT" , ]
  sel <- moy[r,"selection"]
  S2[r,c("surface","PMG","PMG2")] <- S[sel,]
}

R2 <- R2[sort(row.names(R2)),]
S2 <- S2[sort(row.names(S2)),]

H22 <- R2/S2

























# Model lineaire avec planche et NT au hasard -----------------------------

rm(list=ls())


load("../donnees/champ")


# On regarde si pour les non trié, le plot a un effet
{
don <- champ %>% filter(selection == "NT")

mod <- lm(PMG ~ parcelle , data = don)
par(mfrow = c(2,2))
plot(mod)
summary(mod)
anova(mod)

mod <- lm(PMG2 ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)



mod <- lm(surface_recolte_moy ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)


mod <- lm(surface_recolte_moy2 ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)




mod <- lm(surface_recolte_min ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)


mod <- lm(surface_recolte_min2 ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)




mod <- lm(prot_recolte ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)


mod <- lm(poids_max ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)



mod <- lm(poids_min ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)


mod <- lm(GSV ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)

mod <- lm(GSV2 ~ parcelle , data = don)
plot(mod)
summary(mod)
anova(mod)

}

# okay, bon globalement y'a pas d'effet parcelle pour les NT, donc pour les deux NT inconnus, on peut repartir les observations au hasard entre les parcelles pour avoir planche et passage pour toutes les observations et estimer les R




# repartition au hasard entre pas2 pl3 ou pl3 ps2

tmp <- champ %>% filter(planche == "inconnu1" | planche == "inconnu2")

ps2pl3 <- seq(1,nrow(tmp),2)
ps3pl2 <- setdiff(1:nrow(tmp) , ps2pl3)

for (i in ps2pl3){
  tmp[i , c("passage" , "planche")] <- c(2,3)
}

for (i in ps3pl2){
  tmp[i , c("passage" , "planche")] <- c(3,2)
}

don <- champ %>% filter(hauteur > 25)
don[row.names(tmp),] <- tmp

don$selection <- relevel(as.factor(don$selection) , ref = "NT")

graph <- don %>% group_by(planche,passage,selection) %>% summarise() %>% mutate(selection = as.character(selection))

graph$selection2 <- ifelse(graph$selection == "NT" , "t?moin" , graph$selection)

ggplot(graph , aes(x = passage , y = planche , fill = selection)) + geom_tile() + geom_text(aes(label = selection2) , col = "white") + theme(legend.position = "none") + labs(title = "Parcelles au champ")

rm(tmp,graph,i,ps2pl3,ps3pl2)


# quels effets significatifs ?

traits <- c("hauteur" , "nb_epillets" , "surface_recolte_moy" , "surface_recolte_min" , "surface_recolte_max" , "PMG" , "GSV" , "nb_grain" , "surface_recolte_moy2" , "surface_recolte_min2" ,  "surface_recolte_max2" , "PMG2" , "GSV2" , "poids_moy2" , "poids_min2" , "poids_max2" , "prot_recolte")

res <- data.frame()

for (t in traits){
  f <- as.formula(paste(t , "~ selection + passage + planche"))
  
  mod <- lm(f , data = don)
  
  a <- step(mod , direction = "both")
  
  res[t,"effet1"] <- names(a$xlevels)[1]
  res[t,"effet2"] <- names(a$xlevels)[2]
  res[t,"effet3"] <- names(a$xlevels)[3]
  
}

res

# effet de selection : 

traits <- res %>% filter(effet1 == "selection" | effet2 == "selection") %>% row.names()


l <- data.frame()

R <- data.frame()

pval <- data.frame()

for (t in traits){
  ef1 <- res[t , "effet1"]
  ef2 <- res[t , "effet2"]
  
  f <- as.formula(paste(t,"~",ef1,"+",ef2))
  
  mod <- lm(f , data = don)
  a <- drop1(mod , .~. , test = "F")
  
  R[t,"gros"] <- mod$coefficients["selectiongros"]
  R[t,"petit"] <- mod$coefficients["selectionpetit"]
  R[t,"moyen"] <- mod$coefficients["selectionmoyen"]
  
  pval[t,"p-value"] <- a["selection","Pr(>F)"]
  
  lettres <- as.data.frame(cld(glht(model = mod , mcp(selection = "Tukey")) , level = 0.05)$mcletters$Letters)
  
  names(lettres)[1] <- "groups"
  
  lettres$selection <- row.names(lettres)
  
  lettres$variable <- t
  
  lettres$y <- max(don[,t] , na.rm = T) +1
  
  l <- rbind(l , lettres)
  
  # g <- ggplot(don , aes(x=selection , y=don[,t] , col=selection, fill=selection))+
  #   geom_boxplot(alpha = 0.5)+
  #   theme(legend.position="none") +
  #   geom_text(data = lettres, aes(label = groups, y = max(don[,t] , na.rm = T) ), size=5) + stat_mean(col = "black")
  # 
  # print(g)
  # 
  # graphs[[t]] <- g
  
}


# graph et tableau rapport

pander(pval)


don_graph <- don[,c("hauteur","surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","PMG","prot_recolte","selection")] %>% gather(variable,value,1:6)

don_graph$variable2 <- ifelse(don_graph$variable == "hauteur" , "Hauteur" ,
                             ifelse(don_graph$variable == "prot_recolte" , "Taux de protéines",
                                    ifelse(don_graph$variable == "surface_recolte_max2" , "Taille du plus gros grain" , ifelse(don_graph$variable == "surface_recolte_min2" , "Taille du plus petit grain", ifelse(don_graph$variable == "surface_recolte_moy2" , "Taille moyenne des grains" , "PMG")))))

l_graph <- subset(l , variable %in% c("hauteur","surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","PMG","prot_recolte"))

l_graph$variable2 <- ifelse(l_graph$variable == "hauteur" , "Hauteur" ,
                              ifelse(l_graph$variable == "prot_recolte" , "Taux de protéines",
                                     ifelse(l_graph$variable == "surface_recolte_max2" , "Taille du plus gros grain" , ifelse(l_graph$variable == "surface_recolte_min2" , "Taille du plus petit grain", ifelse(l_graph$variable == "surface_recolte_moy2" , "Taille moyenne des grains" , "PMG")))))


don_graph$variable3 <- factor(don_graph$variable2 , levels = c("Taille moyenne des grains" , "Taille du plus petit grain" , "Taille du plus gros grain" , "Hauteur" , "Taux de protéines" , "PMG"))

l_graph$variable3 <- factor(l_graph$variable2 , levels = c("Taille moyenne des grains" , "Taille du plus petit grain" , "Taille du plus gros grain" , "Hauteur" , "Taux de protéines" , "PMG"))

ggplot(don_graph , aes(x=selection , y = value , fill = selection , col = selection)) + geom_boxplot(alpha = 0.5) + stat_mean(col = "black") + facet_wrap(~variable3 , scales = "free") + theme(legend.position = "none") + geom_text(data = l_graph , aes(label = groups, y = y )) + labs(x = "Modalité de sélection" , y = "Valeur des traits" , title = "Effet de la sélection sur différents traits")



tab <- merge(R,pval , by = "row.names") %>% column_to_rownames(var = "Row.names")

tab <- tab[c("surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","hauteur","prot_recolte","PMG"),]

rownames(tab) <- c("Taille moyenne des grains" , "Taille du plus petit grain" , "Taille du plus gros grain" , "Hauteur" , "Taux de protéines" , "PMG")


write.table(round(tab,3) , file = "C:/Users/bienvenu/Documents/Stage_CBienvenu_2023/tables/progres_champ.csv" , sep = ";" , dec = "." , row.names = T , col.names = T)



load("../donnees/S")


R <- as.data.frame(t(R))
S <- S[-4,]

R <- R[sort(row.names(R)) , sort(names(R))]
S <- S[sort(row.names(S)) , sort(names(S))]
R <- R[,-c(1,2,5)]


H2 <- R/S

H2 <- as.data.frame(t(H2))

H2 <- H2[c("surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","PMG"),]

rownames(H2) <- c("Taille moyenne des grains" , "Taille du plus petit grain" , "Taille du plus gros grain" ,"PMG")


H2

write.table(round(H2,3) , file = "C:/Users/bienvenu/Documents/Stage_CBienvenu_2023/tables/H2_champ.csv" , sep = ";" , dec = "." , row.names = T , col.names = T)

ggplot(champ , aes(x=PMG , y=prot_recolte)) + geom_point()
ggplot(champ , aes(x=surface_recolte_moy , y=prot_recolte)) + geom_point()
