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
  
  res[t,"formule"] <- as.character(a$terms)[3]
  
}

res

# effet de selection : 
g <- grep(pattern = "selection" , res$formule)

traits <- traits[g]


l <- data.frame()

R <- data.frame()

pval <- data.frame()

for (t in traits){
  print(t)
  f <- as.formula(paste(t,"~",res[t,"formule"]))
  
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
  
}


# graph et tableau rapport
pval


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



H2

write.table(round(H2,3) , file = "C:/Users/bienvenu/Documents/Stage_CBienvenu_2023/tables/H2_champ.csv" , sep = ";" , dec = "." , row.names = T , col.names = T)

ggplot(champ , aes(x=PMG , y=prot_recolte)) + geom_point()
ggplot(champ , aes(x=surface_recolte_moy , y=prot_recolte)) + geom_point()





















# graph avec tous les R en % de la variance et les p-values et les IC
traits <- c("hauteur" , "nb_epillets" , "surface_recolte_moy2" , "surface_recolte_min2" , "surface_recolte_max" , "PMG" , "GSV" , "nb_grain" , "prot_recolte" )
# on regarde que les sd sont similaires entre les parcelles pour NT
a <- subset(don , selection == "NT")
ec_ty <- data.frame()

for (t in traits){
  f <- as.formula(paste(t,"~ parcelle"))
  b <- aggregate(f , data = a , FUN = sd)
  print(b)
  ec_ty[t,"ecty"] <- mean(b[,t])
}
# ok ça passe


resultat <- data.frame()
i <- 1

for (t in traits){
  
  f <- as.formula(paste(t,"~",res[t,"formule"]))
  
  if (t == "GSV"){f <- as.formula("GSV ~ planche + selection")}
  if (t == "nb_grain"){f <- as.formula("nb_grain ~ planche + selection")}
  
  mod <- lm(f , data = don)
  a <- summary(mod)$coefficients
  
  conf <- confint(mod)
  
  for (sel in c("gros","petit","moyen")){
    chr <- paste0("selection",sel)
    resultat[i,"trait"] <- t
    resultat[i,"R"] <- mod$coefficients[chr] * 100 / ec_ty[t,1]
    resultat[i,"selection"] <- sel 
    resultat[i,"p_value"] <- a[chr,"Pr(>|t|)"]
    resultat[i,"conf_bas"] <- conf[chr,1] * 100 / ec_ty[t,1]
    resultat[i,"conf_haut"] <- conf[chr,2] * 100 / ec_ty[t,1]
    
    i <- i+1
  }

}



resultat$lettres <- ifelse(resultat$p_value < 0.001 , "***" , 
                           ifelse(resultat$p_value < 0.01 , "**",
                                  ifelse(resultat$p_value < 0.05 , "*",NA)))

resultat$trait2 <- ifelse(resultat$trait == "hauteur" , "Hauteur",
                          ifelse(resultat$trait == "nb_epillets" , "Nombre d'épillets",
                                 ifelse(resultat$trait == "nb_grain" , "Nombre de grains par épi",
                                        ifelse(resultat$trait == "prot_recolte" , "Taux N grains",
                                               ifelse(resultat$trait == "surface_recolte_max" , "Taille du plus gros grain",
                                                      ifelse(resultat$trait == "surface_recolte_min2" , "Taille du plus petit grain" , 
                                                             ifelse(resultat$trait == "surface_recolte_moy2" , "Taille moyenne des grains" , resultat$trait)))))))

resultat$trait3 <- factor(resultat$trait2 , levels = c("PMG" , "Taille moyenne des grains" , "Taille du plus petit grain" , "Taille du plus gros grain" , "Hauteur" , "Taux N grains" , "Nombre de grains par épi" , "Nombre d'épillets" , "GSV"))

resultat$htxt <- ifelse(resultat$R < 0 , resultat$conf_bas - 15 , resultat$conf_haut + 5)

ggplot(resultat , aes(x = selection , y = R , fill = selection)) + geom_col() + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3) + geom_text(aes(label = lettres , y = htxt)) + theme(legend.position = "none") + facet_wrap(~trait3) + geom_hline(yintercept = 0) + labs(y = "Progrès estimé (en % de l'écart-type du trait)" , x = "Modalité de sélection" , title = "Progrès estimés après sélection par tamis")


















load("../donnees/opto_recolte_champ")


# repartition au hasard entre pas2 pl3 ou pl3 ps2 et ajout variable selection
for (i in 1:nrow(opto_recolte_champ)){
  p <- opto_recolte_champ$parcelle[i]
  ch <- subset(champ , parcelle == p)
  opto_recolte_champ[i,"selection"] <- unique(ch$selection)
}




row.names(opto_recolte_champ) <- paste(opto_recolte_champ$id , 1:nrow(opto_recolte_champ))




tmp <- opto_recolte_champ %>% filter(planche == "inconnu1" | planche == "inconnu2")

ps <- 2
pl <- 3

tmp[1,"passage"] <- ps
tmp[1,"planche"] <- pl

for (i in 2:nrow(tmp)){
  
  if (tmp[i,"id"] != tmp[i-1,"id"]){
    a <- ps
    ps <- pl
    pl <- a
  }
  
  tmp[i,"passage"] <- ps
  tmp[i,"planche"] <- pl
}


opto_recolte_champ[row.names(tmp),] <- tmp







mod <- lm(Surface ~ selection + passage + planche , data = opto_recolte_champ)

step(mod , direction = "both")

drop1(mod , .~. , test = "F")




summary(mod)
