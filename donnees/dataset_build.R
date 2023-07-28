rm(list = ls())

setwd("~/Stage/donnees")
#setwd("~/Documents/INRA/SelPhen_2023/Stage_CBienvenu_2023/donnees/")

library(tidyverse)


library(nirsextra)
library(rchemo)


# Donnees spectrales ------------------------------------------------------

# Importation des spectres

spectres <- as.data.frame(asd_read_dir("./data_brute/SelPhen_2023_ASD/epo-22-graines"))

p=rbind(list('adj','') , list('red',c(100,100,1))) #, list('detr','')) #,list('snv',''),list('sder',c(1,4,w)))

spectres <- pre(spectres,p)

rm(p)

# separation de l'information du genotype et du numero de grain contenue dans les row.names et creation de deux variables : geno et grain
spectres$geno <- substr(row.names(spectres) , start = 1 , stop = nchar(row.names(spectres))-5)
spectres$grain <- substr(row.names(spectres) , start = nchar(row.names(spectres))-4 , stop = nchar(row.names(spectres)))
spectres$nom_originel <- row.names(spectres)

# homogeneisation a la main des noms de certains genotypes qui sont pas pareils
spectres[which(spectres$geno == "ELAX-35"),"geno"] <- "ELAX_35"
spectres[which(spectres$geno == "ELAX-453"),"geno"] <- "ELAX_453"
spectres[which(spectres$geno == "ELAX-68"),"geno"] <- "ELAX_68"

# tri pour que les numeros de grains aillent dans l'ordre croissant à chaque genotype
spectres <- spectres %>% mutate_at(.vars = "grain" , .funs = as.numeric) %>% relocate(grain , .before = "449") %>% relocate(geno , .before = grain)

# il y a un probleme pour le genotype ELAX_68 : les noms font : 
# ELAX-6800001
# ELAX_6800991
# ELAX_6800992
# ELAX_6800993
# ELAX_6800994
# ELAX_6800995

# ça fait qu'avec la conversion des numeros en dessous, on passe du grain 1 au grain 991 puis 992 etc et c'est pas bien. En plus le ELAX-68 etait un des genotype avec les noms pas homogenes. donc on suppose que le vrai nom de # ELAX-6800001 c'est ELAX_6800990 et on change ça dans les donnees : on remplace grain = 1 par grain = 990

spectres["ELAX-6800001","grain"] <- 990



# il y a un probleme avec le genotype ELAX_35. les noms font :
# ELAX-3500002
# ELAX_3500007
# ELAX_3500008
# ELAX_3500009
# ELAX_3500010
# ELAX_3500011

# La aussi c'etait un des genotypes avec un nom pas homogene. Bon du coup on suppose que le vrai nom de ELAX-3500002 c'est ELAX_3500006 car c'est coherent avec le reste. Donc on remplce grain = 2 par grain = 6 pour que ça fonctionne

spectres["ELAX-3500002","grain"] <- 6


# il y a un probleme avec le genotype ELAX_41 : les noms font :
# ELAX_4101021
# ELAX_4101022
# ELAX_4101024
# ELAX_4101025
# ELAX_4101026

# Donc il en manque 1 (il n'y en a que 5), et on dirait que c'est le ELAX_4101023
# OR, il y a un des genotypes qui s'appelle ELAX_410023 (ce qui donne genotype ELAX_4 grain 10023). A un 1 pres on a ELAX_4101023, le genotype qui manque dans la serie des ELAX_41. Donc on suppose que c'est une erreur de saisie, que le genotype 4 n'existe pas et que c'est bien le 41 manquant. On change ça dans les donnees : 
spectres["ELAX_410023","geno"] <- "ELAX_41"
spectres["ELAX_410023","grain"] <- 1023

# conversion des numeros de grains en 123456 (pour l'instant le numero du grain est donne par le code a 5 chiffre apres le nom du genotype. Par exemple sur ELAX_3500007, le nombre de grain est 7. Comme tout est trie et que ces ce nombre de grain va dans l'ordre croissant pour chaque genotype (sauf les genotypes a probleme la) et bah on retranche a ce nombre de grain la valeur du premier moins 1 comme ça on transforme les valeurs en 1 2 3 4 5 6)

spectres <- spectres %>% arrange(geno,grain)

gen <- spectres[1,"geno"]
n <- spectres[1,"grain"] - 1

for (i in 1:nrow(spectres)){
  if (spectres[i,"geno"] == gen){
    spectres[i , "grain"] <- spectres[i,"grain"] - n}
  
  if (spectres[i,"geno"] != gen){
    gen <- spectres[i,"geno"]
    n <- spectres[i,"grain"] - 1
    spectres[i,"grain"] <- spectres[i,"grain"] - n}
}

summary(spectres$grain)




# On refait tout la meme chose pour spectres bis

spectres_bis <- as.data.frame(asd_read_dir("./data_brute/SelPhen_2023_ASD/R22_epo_BIS_12.07.2022"))

p=rbind(list('adj','') , list('red',c(100,100,1))) #, list('detr','')) #,list('snv',''),list('sder',c(1,4,w)))

spectres_bis <- pre(spectres_bis ,p)

spectres_bis$nom_originel <- row.names(spectres_bis)

# Il y a un probleme sur les GQAX82 et les GQAX4 et les GQAX5
# on suppose qu'ils sont bien dans l'ordre dans le tableau car tous les autres genotypes sont bien dans l'ordre
# modification des GQAX problematiques (il faudra peut-etre les enlever des analyses) :

row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (1)")] <- "GQAX_8200192"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (2)")] <- "GQAX_8200193"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (3)")] <- "GQAX_8200194"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (4)")] <- "GQAX_8200195"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (5)")] <- "GQAX_8200196"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_8200192 (6)")] <- "GQAX_8200197"


row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (10)")] <- "GQAX_400018"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (11)")] <- "GQAX_400019"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (12)")] <- "GQAX_400020"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (7)")] <- "GQAX_400021"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (8)")] <- "GQAX_400022"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_400018 (9)")] <- "GQAX_400023"

row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (1)")] <- "GQAX_500018"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (10)")] <- "GQAX_500019"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (11)")] <- "GQAX_500020"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (12)")] <- "GQAX_500021"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (8)")] <- "GQAX_500022"
row.names(spectres_bis)[which(row.names(spectres_bis) == "GQAX_500018 (9)")] <- "GQAX_500023"


spectres_bis$geno <- substr(row.names(spectres_bis) , start = 1 , stop = nchar(row.names(spectres_bis))-5)
spectres_bis$grain <- substr(row.names(spectres_bis) , start = nchar(row.names(spectres_bis))-4 , stop = nchar(row.names(spectres_bis)))

spectres_bis <- spectres_bis %>% mutate_at(.vars = "grain" , .funs = as.numeric) %>% relocate(grain , .before = "449") %>% relocate(geno , .before = grain) %>% arrange(geno,grain)
# pas de probleme de noms à homogeneiser ici

# pour spectres bis on numerote les genotypes de 7 a 12 au lieu de 1bis à 6bis, ce sera plus pratique pour la suite
gen <- spectres_bis[1,"geno"]
n <- spectres_bis[1,"grain"] - 7

for (i in 1:nrow(spectres_bis)){
  if (spectres_bis[i,"geno"] == gen){
    spectres_bis[i , "grain"] <- spectres_bis[i,"grain"] - n}
  
  if (spectres_bis[i,"geno"] != gen){
    gen <- spectres_bis[i,"geno"]
    n <- spectres_bis[i,"grain"] - 7
    spectres_bis[i,"grain"] <- spectres_bis[i,"grain"] - n}
}

summary(spectres_bis$grain)



# Assemblage final : on reunit tout dans un tableau
spectres <- rbind(spectres , spectres_bis) %>% relocate(nom_originel , .after = grain)
summary(spectres$grain)
# tout ça m'a l'air plutôt bon

# On enleve l'information inutile et redondante du nom des genotypes car ce sera plus simple pour faire correspondre avec les donnees de l'optomachine

for (i in 1:nrow(spectres)){
  if (!spectres[i,"geno"] %in% c("ANVERGUR" , "IXOS" , "OBELIX") ){
    spectres[i,"geno"] <- strsplit(spectres[i,"geno"] , split = "_")[[1]][2]
  }
}


# on trie et renomme pour que ce soit clair
spectres <- spectres %>% arrange(geno,grain)

row.names(spectres) <- paste0(spectres$geno,"_",spectres$grain)

names(spectres)[4:ncol(spectres)] <- paste0("X" , names(spectres)[4:ncol(spectres)])

#sauvegarde du tableau obtenu (avec les colonnes geno et grain  et nom originel au cas où)
spectres_info <- spectres
save(spectres_info , file = "spectres_info")


# Sauvegarde d'un autre tableau avec seulement les longueurs d'onde pour que ce soit plus simple à utiliser

spectres_brutes <- spectres[-c(1:3)]

# representation des spectres

#plotsp(spectres_brutes)

# dejà y'a un spectre dissident, on l'enleve

sp_brutes <- subset(spectres_brutes , X1000 > 0.1)

save(sp_brutes , file = "sp_brutes")


# Prétraitements

load("sp_brutes")

sp_SG <- savgol(sp_brutes , m = 0 , n = 15 , p = 1)

sp_det <- as.data.frame(detrend(sp_SG , degree = 1))

sp_det_SNV <- as.data.frame(snv(sp_det))

sp_dev1 <- as.data.frame(dderiv(sp_SG , n = 15))

sp_dev2 <- as.data.frame(dderiv(sp_dev1 , n = 21))


# sauvegarde des spectres traites
spectres <- list(sp_brutes , sp_det , sp_det_SNV , sp_dev1 ,sp_dev2)
names(spectres) <- c("brutes","det","det_SNV","dev1","dev2")

save(spectres , file = "spectres")

rm(list = ls())





# On fait aussi une liste avec les spectres moyens de chaque génotype

load("spectres")

sp <- spectres$brutes
sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
brutes <- sp %>% group_by(geno) %>% summarise_all(.funs = mean) %>% column_to_rownames(var = "geno")

sp <- spectres$det
sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
det <- sp %>% group_by(geno) %>% summarise_all(.funs = mean) %>% column_to_rownames(var = "geno")

sp <- spectres$det_SNV
sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
det_SNV <- sp %>% group_by(geno) %>% summarise_all(.funs = mean) %>% column_to_rownames(var = "geno")

sp <- spectres$dev1
sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
dev1 <- sp %>% group_by(geno) %>% summarise_all(.funs = mean) %>% column_to_rownames(var = "geno")

sp <- spectres$dev2
sp$geno <- sapply(strsplit(row.names(sp) , "_") , "[" , 1)
dev2 <- sp %>% group_by(geno) %>% summarise_all(.funs = mean) %>% column_to_rownames(var = "geno")


spectres_moy <- list(brutes , det , det_SNV , dev1 , dev2)
names(spectres_moy) <- c("brutes","det","det_SNV","dev1","dev2")

save(spectres_moy , file = "spectres_moyens")

rm(list = ls())


# Donnees optomachine grains semes -----------------------------------------------------

rm(list=ls())

# recuperation de tous les noms
file_names <- list.files(path = "./data_brute/SelPhen_2023_Optoagri")

# Initialisation du tableau
opto <- read.table(file = "./data_brute/SelPhen_2023_Optoagri/2022-11-02_141909_0101_Mecarnot_JRL_1/2022-11-02_141909_0101_Mecarnot_JRL_1_METRO.dat" , header = T , sep = "\t" , dec = ",",fileEncoding="latin1")

opto$geno <- "JRL_1"


# Remplissage du tableau
for (i in file_names[2:61]){
  df <- read.table(file = paste0("./data_brute/SelPhen_2023_Optoagri/" , i , "/" , i , "_METRO.dat") , header = T , sep = "\t" , dec = ",",fileEncoding="latin1")
  
  df$geno <- strsplit(i , split = "t_")[[1]][2]
  
  opto <- rbind(opto , df)
}

opto <- opto %>% relocate(geno , .before = Réf..Ech)

# formattage du tableau pour les genotypes marques de façon pas homogene (construction au cas par cas)

# donnees JRL_1
g<- c(3,4,5,6,9,10)
for (i in which(opto$geno == "JRL_1")){
  if (opto[i,"Index"] <= 6){
    opto[i,"geno"] <- g[1]
    opto[i,"grain"] <- opto[i,"Index"]}
  
  if (opto[i,"Index"] > 6 & opto[i,"Index"] <= 12){
    opto[i,"geno"] <- g[2]
    opto[i,"grain"] <- opto[i,"Index"] - 6}
  
  if (opto[i,"Index"] > 12 & opto[i,"Index"] <= 18){
    opto[i,"geno"] <- g[3]
    opto[i,"grain"] <- opto[i,"Index"] - 12}
  
  if (opto[i,"Index"] > 18 & opto[i,"Index"] <= 24){
    opto[i,"geno"] <- g[4]
    opto[i,"grain"] <- opto[i,"Index"] - 18}
  
  if (opto[i,"Index"] > 24 & opto[i,"Index"] <= 30){
    opto[i,"geno"] <- g[5]
    opto[i,"grain"] <- opto[i,"Index"] - 24}
  
  if (opto[i,"Index"] > 30){
    opto[i,"geno"] <- g[6]
    opto[i,"grain"] <- opto[i,"Index"] - 30}
} 

opto <- opto %>% relocate(grain , .before = Réf..Ech)

# donnees 16-70
g<- c(16,27,28,49,51,70)
for (i in which(opto$geno == "16-70")){
  if (opto[i,"Index"] <= 6){
    opto[i,"geno"] <- g[1]
    opto[i,"grain"] <- opto[i,"Index"]}
  
  if (opto[i,"Index"] > 6 & opto[i,"Index"] <= 12){
    opto[i,"geno"] <- g[2]
    opto[i,"grain"] <- opto[i,"Index"] - 6}
  
  if (opto[i,"Index"] > 12 & opto[i,"Index"] <= 18){
    opto[i,"geno"] <- g[3]
    opto[i,"grain"] <- opto[i,"Index"] - 12}
  
  if (opto[i,"Index"] > 18 & opto[i,"Index"] <= 24){
    opto[i,"geno"] <- g[4]
    opto[i,"grain"] <- opto[i,"Index"] - 18}
  
  if (opto[i,"Index"] > 24 & opto[i,"Index"] <= 30){
    opto[i,"geno"] <- g[5]
    opto[i,"grain"] <- opto[i,"Index"] - 24}
  
  if (opto[i,"Index"] > 30){
    opto[i,"geno"] <- g[6]
    opto[i,"grain"] <- opto[i,"Index"] - 30}
} 


# donnees 68-59
g<- c(68,69,67,65,41,59)
for (i in which(opto$geno == "68-59")){
  if (opto[i,"Index"] <= 6){
    opto[i,"geno"] <- g[1]
    opto[i,"grain"] <- opto[i,"Index"]}
  
  if (opto[i,"Index"] > 6 & opto[i,"Index"] <= 12){
    opto[i,"geno"] <- g[2]
    opto[i,"grain"] <- opto[i,"Index"] - 6}
  
  if (opto[i,"Index"] > 12 & opto[i,"Index"] <= 18){
    opto[i,"geno"] <- g[3]
    opto[i,"grain"] <- opto[i,"Index"] - 12}
  
  if (opto[i,"Index"] > 18 & opto[i,"Index"] <= 24){
    opto[i,"geno"] <- g[4]
    opto[i,"grain"] <- opto[i,"Index"] - 18}
  
  if (opto[i,"Index"] > 24 & opto[i,"Index"] <= 30){
    opto[i,"geno"] <- g[5]
    opto[i,"grain"] <- opto[i,"Index"] - 24}
  
  if (opto[i,"Index"] > 30){
    opto[i,"geno"] <- g[6]
    opto[i,"grain"] <- opto[i,"Index"] - 30}
} 


# donnees 31-13
g<- c(31,29,38,40,35,13)
for (i in which(opto$geno == "31-13")){
  if (opto[i,"Index"] <= 6){
    opto[i,"geno"] <- g[1]
    opto[i,"grain"] <- opto[i,"Index"]}
  
  if (opto[i,"Index"] > 6 & opto[i,"Index"] <= 12){
    opto[i,"geno"] <- g[2]
    opto[i,"grain"] <- opto[i,"Index"] - 6}
  
  if (opto[i,"Index"] > 12 & opto[i,"Index"] <= 18){
    opto[i,"geno"] <- g[3]
    opto[i,"grain"] <- opto[i,"Index"] - 12}
  
  if (opto[i,"Index"] > 18 & opto[i,"Index"] <= 24){
    opto[i,"geno"] <- g[4]
    opto[i,"grain"] <- opto[i,"Index"] - 18}
  
  if (opto[i,"Index"] > 24 & opto[i,"Index"] <= 30){
    opto[i,"geno"] <- g[5]
    opto[i,"grain"] <- opto[i,"Index"] - 24}
  
  if (opto[i,"Index"] > 30){
    opto[i,"geno"] <- g[6]
    opto[i,"grain"] <- opto[i,"Index"] - 30}
} 


# Pour les donnees 3_10_bis, 13_35_bis, 38_51_bis, et 59_70_bis on ne peut pas savoir qui est qui donc on ne prend pas ces donnees

opto <- subset(opto , !geno %in% c("3_10_bis" , "16_35_bis" , "38_51_bis" , "59_70_bis"))

row.names(opto) <- 1:nrow(opto)

# Construction du tableau pour les genotypes marques de façon homogene

for (i in 145:2014){
  if (opto[i,"Index"] %in% 1:12){
    opto[i,"geno"] <- strsplit(opto[i,"geno"] , split = "_")[[1]][1]
    opto[i,"grain"] <- opto[i,"Index"]
  }
  
  if (opto[i,"Index"] %in% 13:24){
    opto[i,"geno"] <- strsplit(opto[i,"geno"] , split = "_")[[1]][2]
    opto[i,"grain"] <- opto[i,"Index"] - 12
  }
  
  if (opto[i,"Index"] %in% 25:36){
    opto[i,"geno"] <- strsplit(opto[i,"geno"] , split = "_")[[1]][3]
    opto[i,"grain"] <- opto[i,"Index"] - 24
  }
}

summary(opto$grain)


# Pour le genotype 453 : 

opto[which(opto$geno == 453),"grain"] <- opto[which(opto$geno == 453),"Index"]

# Changement des noms de variables pour que ce soit plus clair
names(opto) <- c("geno",                      "grain"             ,        "Ref.Ech"  ,              
                 "Index",                     "Longueur"      ,       "Longueur.interieure" ,
                 "Largeur",              "Perimetre"      ,      "Perimetre.de.Crofton",
                 "Perimetre.convexe",    "Dimetre.Eq"     ,     "Surface"            ,
                 "Surface.convexe"  ,   "Surface.DiffCAP" ,    "Surface.DiffCEN"    ,
                 "Surface.DiffEAP"  ,   "Finesse"                ,   "Excentricite"             ,
                 "F.Feret"                ,  "Compacite"               ,  "Circularite"              ,
                 "Rugosite"                 , "Index.de.courbure",         "Moment.inertie.1"         ,
                 "Moment.inertie.2",          "Moment.inertie.3"  ,        "Symetrie"                 ,
                 "Moyenne.ndg"      ,         "Ecart.type.ndg"     ,       "Minimum.ndg"              ,
                 "Maximum.ndg"       ,        "Histo.Kurtose"      ,      "Histo.Moyenne"           ,
                 "Histo.Pic"         ,       "Histo.Assymetrie"    ,     "Histo.Ecart.type"        ,
                 "Histo.Variance"     ,      "Cooc.Uniformite"      ,    "Cooc.Contraste"          ,
                 "Cooc.Correlation"    ,     "Cooc.Variance.globale" ,   "Cooc.Homogeneïte"        ,
                 "Cooc.Somme.des.moyennes",  "Cooc.Somme.des.variances", "Cooc.Somme.des.entropies",
                 "Cooc.Entropie.globale",    "Cooc.Ecart.variance",      "Cooc.Ecart.entropie"     ,
                 "Cooc.Correlation.IC1"  ,   "Cooc.Correlation.IC2",     "RVB.R.Moy"                ,
                 "RVB.R.Ect"               ,  "RVB.R.Min"             ,    "RVB.R.Max"                ,
                 "RVB.V.Moy",                 "RVB.V.Ect"              ,   "RVB.V.Min"                ,
                 "RVB.V.Max" ,                "RVB.B.Moy"               ,  "RVB.B.Ect"                ,
                 "RVB.B.Min"  ,               "RVB.B.Max"                , "TSI.T.Moy"                ,
                 "TSI.T.Ect"   ,              "TSI.T.Min",                 "TSI.T.Max"                ,
                 "TSI.S.Moy"    ,             "TSI.S.Ect" ,                "TSI.S.Min"                ,
                 "TSI.S.Max"     ,            "TSI.I.Moy"  ,               "TSI.I.Ect"                ,
                 "TSI.I.Min"      ,           "TSI.I.Max"   ,              "CMJ.C.Moy"                ,
                 "CMJ.C.Ect"       ,          "CMJ.C.Min"    ,             "CMJ.C.Max"                ,
                 "CMJ.M.Moy"        ,         "CMJ.M.Ect"     ,            "CMJ.M.Min"                ,
                 "CMJ.M.Max"         ,        "CMJ.J.Moy"      ,           "CMJ.J.Ect"                ,
                 "CMJ.J.Min"          ,       "CMJ.J.Max"       ,          "Lab.L.Moy"                ,
                 "Lab.L.Ect"           ,      "Lab.L.Min"        ,         "Lab.L.Max"                ,
                 "Lab.a.Moy"            ,     "Lab.a.Ect"         ,        "Lab.a.Min"                ,
                 "Lab.a.Max"             ,    "Lab.b.Moy"          ,       "Lab.b.Ect"                ,
                 "Lab.b.Min"              ,   "Lab.b.Max"           ,      "X" )


# On met les row.names qui vont bien
row.names(opto) <- paste0(opto$geno,"_",opto$grain)

# les colonne ref ech, X et index servent à rien, on les enlève
unique(opto$Ref.Ech)
unique(opto$X)

opto <- opto %>% select(!c("Ref.Ech" , "X" , "Index"))

# Sauvegarde du tableau
save(opto , file = "opto")

rm(g,i,df,file_names,opto)


# Donnees des bacs --------------------------------------------------------
{
  
  rm(list = ls())
  
  
  # Chargement des donnees des bacs :
  bac <- read.table("./data_brute/BAC1.csv" , header = T , sep = ";" , dec = ".")
  
 
  
  for (i in 2:6){
    tmp <- read.table(paste0("./data_brute/BAC",i,".csv") , header = T , sep = ";" , dec = ".")
    names(tmp) <- names(bac)
    bac <- rbind(bac , tmp)
  }
  
  rm(tmp , i)
  
  # Modification de la variable grain en fonction de l'information dans la variable graines pour que ce soit homogene avec le reste (de 1 a 12)
  
  bac$grain <- ifelse(bac$graines == "01:06" , bac$grain , bac$grain + 6)
  
  
  # creation d'une variable geno homogene aux autres tableaux de donnees
  # il y a des noms avec .sp à la fin, on les enleve, et on separe le ELAX du numero
  
  a <- strsplit(bac$gen , split = ".sp")
  l <- sapply(a, "[",1)
  b <- strsplit(l , split = "_")
  bac$geno <- sapply(b, "[",2)
  
  as.numeric(bac$geno) # c'est bon on a que des numeros (pas de message d'erreur)
  
  rm(a,b,l)
  
  
  # On cree le row.names pour que ce soit homogene aux autres tableaux
  
  row.names(bac) <- paste0(bac$geno , "_" , bac$grain)
  # comme les noms doivent etre unique, on voit qu'il n'y a pas de grain qui apparaisse 2 fois dans le tableau donc c'est bon signe
  
  
  # modifications pour les genotypes qui ont pas germe
  # On commence par creer une variable qui indique l'individu semé au premier semis
  
  bac$semis_1 <- row.names(bac)
  
  # Les individus qui ont pas germe sont ressences dans le fichier "correction_bac.csv"
  # importation de ce fichier
  correc <- read.table("./data_brute/correction_bac.csv" , header = T , sep = ";" , dec = ".")
  
  
  
  # y'a une poignee d'individus pour lesquels l'info dans les donnees numerisee est pas bonne. ce sont les individus surlignes en rose du bac 2 dans la version papier.
  # Donc on fait les corrections à la main.
  # On modifie aussi les noms dans correc pour que ceux ci soient quand meme traites s'ils ont pas germe
  
  # pour ceux où il y a une erreur dans la saisie et qui n'ont pas germe
  for (i in c("311" , "222" , "212" , "69")){
    bac[paste0(i,"_2") , "grain"] <- 12
    row.names(bac)[which(row.names(bac)==paste0(i,"_2"))] <- paste0(i,"_12")
    correc[which(correc$ind == paste0(i,"_2")) , "ind"] <- paste0(i,"_12")
  }
  
  
  # pour ceux où il y a juste une erreur dans la saisie mais qui ont germe on a pas besoin de modifier correc
  
  for (i in c("489" , "124" , "316" , "341" , "233" , "399" , "166" , "240")){
    bac[paste0(i,"_2") , "grain"] <- 12
    row.names(bac)[which(row.names(bac)==paste0(i,"_2"))] <- paste0(i,"_12")
  }
  
  # le fichier correc ressence tous les individus qui n'ont pas germe, et la colonne de la boite dans laquelle le grain de remplacement a ete pris. Dans ce tableau il y a tous les individus pour lesquels le grain de remplacement a ete pris dans la boîte bis correspondant à la boite de base (la boite de base qui contenait la graine qui n'a pas germe, la boite bis qui contenait la graine qui a remplace celle qui n'a pas germe) c'est à dire meme numero de boîte et meme rang (donc meme genotype). Comme le numero de boite et le rangs definissent le genotype, on a besoin seulement de la colonne pour savoir quel grain a ete pris (car le numero de boite et le rang correspondent entre les boites de base et les boites bis) : 
  # Du coup dans la boîte de base la disposition est telle que les grains dans les colonnes 1 à 6 sont bien les grains 1 à 6, et les grains dans les colonnes 7 à 12 sont aussi des grains 1 à 6.
  # dans la boîte bis, la disposition est telle que les grains dans les colonnes 7 à 12 sont bien les grains 7 à 12, et les grains dans les colonnes 1 à 6 sont des grains 7 à 12.
  # Donc pour savoir le numero de grain : si c'est entre 7 et 12 c'est bon, si c'est entre 1 et 6, il faut rajouter +6 pour avoir le bon numero de grain. (car depuis le debut on a definit que 1 à 6 ça reste 1 à 6, mais de 1bis à 6bis, ça devient 7 à 12).
  # pour certains genotypes on sait pas par quoi ils ont ete remplaces donc ils ont un point d'interrogtion
  
  correc$col <- as.numeric(correc$col) # points d'interrogation transformes en NA
  correc$grain <- ifelse(correc$col >= 7 , correc$col , correc$col + 6)
  
  summary(correc$grain) # c'est good
  
  # on cree une colonne new_ind qui donne l'identifiant de l'individu qui a ete plante en remplacement
  
  correc$geno <- sapply(strsplit(correc$ind , split = "_") , "[",1)
  correc$new_ind <- paste0(correc$geno , "_" , correc$grain)
  
  row.names(correc) <- correc$ind
  
  
  # on modifie bac : on crée une colonne semis 2 qui va ressencer les individus semés au deuxième semis
  c <- 1 # compteur utile dans le if des genotypes inconnus
  bac$semis_2 <- NA
  
  for (i in 1:nrow(correc)){
    ind <- correc[i,"ind"]
    new_ind <- correc[i,"new_ind"]
    
    # Si on sait quel genotype a remplace : on remplit la colonne semis_2
    if (is.na(correc[ind , "grain"]) == FALSE){
      bac[ind,"semis_2"] <- new_ind
    }
    
    # Si on ne sait pas quel genotype a ete seme (les NA dans correc), on met INCONNU dans semis_2
    
    if (is.na(correc[ind , "grain"]) == TRUE){
      bac[ind,"semis_2"] <- paste0("INCONNU_",c)
      c <- c+1
    }
  }

  rm(ind , new_ind , i)
  
  
  
  # ça c'etait pour toutes les donnees surlignees en orange sur le papier sans annotations roses (sauf pour la boite 12 où c'est pris en compte malgre les annotations roses).
  
  # il y a d'autres donnees plus chiantes qu'on peut pas simplement remplacer comme ça car les boites et les rangs entre l'individu remplace et le remplaçant ne correspondent pas. Du coup faut aller chercher à la main à quels genotype ils correspondent. Ce sont toutes les donnees entourrees en rose sur le format papier (sauf celles de la boite 12)
  
  
  func <- function(b,r,d,g,ind){
    # recuperation des infos sur le genotype remplaçant
    a <- subset(bac , boite == b & Rang.OK == r & Debut == d & grain == g)
    
    # genotype remplace = ind
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <<- bac[ind , ]
    bac[new_ind,"geno"] <<- new_geno
    bac[new_ind,"grain"] <<- new_grain
    bac[new_ind,"semis_1"] <<- ind
    bac[new_ind,"semis_2"] <<- new_ind
    bac <<- bac[-which(bac$semis_1==ind & is.na(bac$semis_2)==T),]
  }
  
  #func(b=4 , r="E" , d=1 , g=1 , ind="394_1")
  func(b=6 , r="E" , d=3 , g=3 , ind="305_2")
  func(b=10 , r="F" , d=3 , g=3 , ind="488_7")
  func(b=10 , r="C" , d=10 , g=4 , ind="188_2")
  func(b=2 , r="C" , d=8 , g=2 , ind="65_2")
  func(b=10 , r="E" , d=8 , g=2 , ind = "447_3")
  #func(b=9 , r="D" , d=7 , g=1 , ind="428_5")
  func(b=8 , r="C" , d=2 , g=2 , ind="428_6")
  func(b=7 , r="H" , d=2 , g=2 , ind="311_6")
  func(b=4 , r="A" , d=1 , g=1 , ind="489_6")
  #func(b=1 , r="B" , d=8 , g=2 , ind="6_6")
  #func(b=1 , r="B" , d=8 , g=2 , ind="6_7")
  
  #pb pour le 428_5, le 394_1, le 6_6 et le 6_7 : genotype de rempplacement deja plante
  # on les remplace par des genotypes inconnus
 
  
  for (g in c("428_5" , "394_1" , "6_6" , "6_7")){
    bac[g,"semis_2"] <- paste0("INCONNU_",c)
    bac[g,"geno"] <- "INCONNU"
    bac[g,"grain"] <- c
    c <- c+1
  }
  
  
  # On actualise les rownames
  row.names(bac) <- ifelse(is.na(bac$semis_2)==T,bac$semis_1,bac$semis_2)
  
  # On enleve les variables inutiles 
  bac$X.1 <- bac$Fin <- bac$boite <- bac$boite.1 <- bac$Debut <- bac$RANG <- bac$Rang.OK <- bac$graines <- bac$gen <- bac$grain <-NULL
  
  # on crée la variable de date de semis utile
  bac$semis <- ifelse(is.na(bac$semis_2) == T , "06/01" , "17/02")

  save(bac , file = "bac")
  
  rm(bac,correc,c,g,func)
  
}



# Ajout des donnees d'epiaison
{
load("bac")
epiaison <- read.table("data_brute/epiaison.csv" , sep = ";" , dec = "." , header = T , na.strings = "")

bac$key <- paste0(bac$BAC,"_",bac$X,"_",bac$Y)

epiaison$key <- paste0(epiaison$BAC,"_",epiaison$X,"_",epiaison$Y)
epiaison <- epiaison %>% select(key,DATE)

bac <- merge(bac , epiaison , by = "key") %>% select(!c("key")) %>% rename(epiaison = DATE)

row.names(bac) <- ifelse(is.na(bac$semis_2)==T,bac$semis_1,bac$semis_2)


# Les donn?es d'?piaison sont cod?es comme le jour de mai (2 = 2 mai)
# Le premier semis a ?t? fait le 6 janvier, et le deuxi?me semis a ?t? fait le 17 f?vrier
# On corrige la colonne d'epiaison en mettant le nombre de jours r?els avant epiaison
# pour le semis 1, il y a 115 jours entre le 6 janvier et le 1 mai, donc on rajoute 115 au donnees du premier semis
# pour le semis 2, il y a 63 jors entre le 17 fevrier et le 1 mai, donc on rajoute 60 aux donnees du deuxieme semis

bac$epiaison <- as.numeric(bac$epiaison)

bac$epiaison <- ifelse(is.na(bac$semis_2)==T , bac$epiaison + 115 , bac$epiaison + 63)

# On transforme en degre jour
meteo <- read.table(file = "data_brute/essai_INRAE_STATION_34172005 (14).csv" , sep = ";" , dec = "." , header = T , skip = 9)

# degre jours depuis le 6 janvier
dj <- meteo %>% filter(MOIS > 1 | JOUR >= 6)
row.names(dj) <- 1:nrow(dj)
dj$DJ <- NA
dj$DJ[1] <- dj$TM[1]

for (i in 2:nrow(dj)){
  dj[i,"DJ"] <- dj[i-1,"DJ"] + dj[i,"TM"]
}


# degre jour depuis le 17 fevrier
a <- dj[which(dj$MOIS == 2 & dj$JOUR == 16) , "DJ"]

dj$DJ2 <- dj$DJ - a

dj2 <- dj %>% filter(DJ2 >0)


# On inclue ça dans la variable epiaison
bac$preco <- NA

for (i in 1:nrow(bac)){
  index <- bac[i,"epiaison"]
  if (bac[i,"semis"] == "06/01"){
    bac[i,"preco"] <- dj[index,"DJ"]
  }
  
  if (bac[i,"semis"] == "17/02"){
    bac[i,"preco"] <- dj2[index,"DJ2"]
  }
}
}

# Ajout des donnees de luzerne
{
bac$luz <- ifelse(bac$BAC==1 | bac$BAC==4 , "catera" ,
                  ifelse(bac$BAC==6 | bac$BAC==2 , "aria" , "sans"))


save(bac , file = "bac")

rm(list=ls())
}


# Ajout des données d'azote de feuille drapeau
{
load("bac")
N <- read.table("./data_brute/SelPhen2023_azote_flo_drapeau_azmod4_3.txt" , header = T , sep = "," , dec = ".")

# formatage du tableau
N$BAC2 <- sapply(str_split(N$list , "_") , "[" , 1)
N$tmp <- sapply(str_split(N$list , "_") , "[" , 2)
N$X <- as.numeric(sapply(str_split(N$tmp , "000") , "[" , 2))
N$Ytmp <- sapply(str_split(N$tmp , "000") , "[" , 1)
N$Y <- as.numeric(sapply(str_split(N$Ytmp , "y") , "[" , 2))

N$list <- N$tmp <- N$Ytmp <- NULL


# Il y a un problem avec les donnes les lignes y10, elles sont notees 1. On resout ce probleme

for (i in 2:nrow(N)){
  if ((N$Y[i-1] == 9 & N$Y[i] == 1) | (N$Y[i-1] == 10 & N$Y[i] == 1) ){
    N$Y[i] <- 10
  }
}


# On trouve pour quels bacs il manque des données
N$count <- 1
cherche <- N %>% group_by(BAC2 , Y) %>% summarise(count=sum(count)) %>% filter(count != 32)


# On regarde a la main pour essaye de trouver quelles donnees il manque exactement en prenant les plantes manquantes comme repere
# 
# alamain <- N %>% filter((BAC2 == cherche$BAC2[1] & Y == cherche$Y[1]) |
#                           (BAC2 == cherche$BAC2[2] & Y == cherche$Y[2]) |
#                           (BAC2 == cherche$BAC2[3] & Y == cherche$Y[3]) |
#                           (BAC2 == cherche$BAC2[4] & Y == cherche$Y[4]) |
#                           (BAC2 == cherche$BAC2[5] & Y == cherche$Y[5]) |
#                           (BAC2 == cherche$BAC2[6] & Y == cherche$Y[6]) |
#                           (BAC2 == cherche$BAC2[7] & Y == cherche$Y[7]) |
#                           (BAC2 == cherche$BAC2[8] & Y == cherche$Y[8]))
# 
# 
# 
# # cellule de debug
# i <- 208
# x <- 1
# {
#   print(alamain[c(i,i+1),])
#   i <- i+2
#   x <- x+1
# }

# Alright

# bac x9y5 : y=10 il manque x=16
# bac x9y5 : y=13 il manque x=1 et x=2
# bac x9y4 : y=4 il manque x=31 et x=32
# bac x9y4 : y=5 il manque x=1 et x=2
# bac x9y4 : y=9 il manque x=31 et x=32
# bac x9y4 : y=13 il manque x=1 x=2 x=3 x=4 x=31 x=32
# bac x10y5 : y=3 il manque x=1 x=2
# bac x10y5 : y=11 il manque x=1 x=2


# Mais du coup c'est bon y'a rien a faire care les manquants sont pas dans les donnees et ce qui ont ete mesures ont bien les bons x et y donc c'est bon le merge va tou bien faire



# Les donnees ont ete acquises lignes y par lignes y en faisant 32 mesures : 2 par plante; Donc il faut transformer la variable X pour qu elle devienne vraie

N$Xfaux <- N$X

N$X <- N$Xfaux - N$Xfaux %/% 2


# On fait la moyenne des deux observations pour chaque individu

N2 <- N %>% group_by(BAC2 , Y , X) %>% summarise(N_flag = mean(az))


# Transposition de la variable bac

N2$BAC <- ifelse(N2$BAC2 == "x09y03" , "1" ,
                ifelse(N2$BAC2 == "x09y04" , "2" ,
                       ifelse(N2$BAC2 == "x09y05" , "3" ,
                              ifelse(N2$BAC2 == "x10y03" , "4" ,
                                     ifelse(N2$BAC2 == "x10y04" , "5" , "6")))))




# On enlève les donnees qui ont ete prises dans l'air

la <- read.table("./data_brute/feuille_ou_air.txt" , header = T , sep = ",")

la$BAC2 <- sapply(strsplit(la$list , "_") , "[" , 1)
la$pos <- sapply(strsplit(la$list , "_") , "[" , 2)
la$pos2 <- sapply(strsplit(la$pos , "y") , "[" , 2)
for (i in 1:nrow(la)){
  la[i,"Y"] <- substr(la[i,"pos2"] , start = 1 , stop = nchar(la[i,"pos2"])-5)
  
  la[i,"Xfaux"] <- as.numeric(substr(la[i,"pos2"] , start = 4 , stop = nchar(la[i,"pos2"])))
}

la$X <- la$Xfaux - la$Xfaux %/% 2

la2 <- la %>% group_by(BAC2 , Y , X) %>% summarise(id = mean(id))

table(la2$id)

la2[which(la2$id == 0.5),"id"] <- 0

table(la2$id)

la2$BAC <- ifelse(la2$BAC2 == "x09y03" , "1" ,
                 ifelse(la2$BAC2 == "x09y04" , "2" ,
                        ifelse(la2$BAC2 == "x09y05" , "3" ,
                               ifelse(la2$BAC2 == "x10y03" , "4" ,
                                      ifelse(la2$BAC2 == "x10y04" , "5" , "6")))))



# Ajout dans les donnees bac

la2$key <- paste0(la2$BAC , "_" , la2$X , "_" , la2$Y)
bac$key <- paste0(bac$BAC , "_" , bac$X , "_" , bac$Y)
N2$key <- paste0(N2$BAC , "_" , N2$X , "_" , N2$Y)

N2$X <- N2$Y <- N2$BAC <- la2$X <- la2$Y <- la2$BAC <- la2$BAC2 <- NULL

bac <- merge(bac , N2 , by = "key" , all = T)

bac <- merge(bac , la2 , by = "key" , all = T)

bac$key <- NULL

row.names(bac) <- ifelse(is.na(bac$semis_2)==T,bac$semis_1,bac$semis_2)

bac$N_flag <- ifelse(bac$id == 0 , NA , bac$N_flag)

bac$id <- NULL


save(bac , file = "bac")

rm(list=ls())
}



# ajout des donnees de teneur en proteine des grains plantes
{
load("bac")

pred1 <- read.table(file = "data_brute/Pred_prot_grain_ASD_nov2022.csv" , sep = "," , dec = "." , header = T)

pred2 <- read.table(file = "data_brute/Pred_prot_grain_ASD_nov2022_bis.csv" , sep = "," , dec = "." , header = T)

pred <- rbind(pred1,pred2) %>% mutate(X = NULL)

rm(pred1,pred2)

# On ressort le jeu de donnees spectres_info qui contient le lien entre les noms et l'identifiant du grain
load("spectres_info")

spectres_info <- spectres_info %>% select(nom_originel , grain , geno) %>% rename(nom = nom_originel)


# verification qu'on a bien la même info :

unique(sort(pred$nom)==sort(spectres_info$nom)) 
# c'est good


tout <- merge(pred,spectres_info , by = "nom")

row.names(tout) <- paste0(tout$geno,"_",tout$grain)


final <- tout %>% rename(prot_semis = prot_calib_grain) %>% select(prot_semis)


bac <- merge(bac , final , by = "row.names" , all.x = T) %>% column_to_rownames(var = "Row.names")

bac$BAC <- as.factor(bac$BAC)

save(bac , file = "bac")


load("opto")

opto <- merge(opto , final , by = "row.names" , all.x = T) %>% column_to_rownames(var = "Row.names")

save(opto , file = "opto")

rm(list = ls())
}




## Ajout des donnees de taille du grain

load("bac")

load("opto")

a <- opto %>% select("Surface")

bac <- merge(bac,a,by = "row.names" , all.x = T) %>% column_to_rownames(var = "Row.names")

save(bac , file = "bac")



# Ajout des donnees de hauteur du brin maitre

rm(list = ls())

load("bac")

h <- read.table("./data_brute/mesure_hauteur.csv" , header = F , sep = ";" , dec = ".")
names(h) <- c("geno" , "hauteur")

h <- h %>% column_to_rownames(var = "geno")
# mince y'a des genotypes sur lesquels on s'est trompe. On les enlève

h <- h[-which(h$geno == "112_5" | h$geno == "130_5" | h$geno == "INCONNU_10"),]


row.names(h) <- h$geno
h$geno <- NULL

bac <- merge(bac , h , by = "row.names" , all.x = T) %>% column_to_rownames(var = "Row.names")

bac$hauteur <- as.numeric(bac$hauteur)
bac$hauteur <- bac$hauteur + 30


save(bac , file = "bac")





# ajout des donnees poids d'epi

rm(list = ls())

load("bac")


p <- read.table("./data_brute/POIDS_EPI_CLEM.csv" , header = F , sep = ";")


# y'a que 1 colonne où les lignes vont 3 par 3, 1:nom du genotype, 2:nb epi, 3:poids epi en mg

tab <- data.frame()

index <- seq(1,nrow(p),3)

for (i in index){
  tab[p[i,1],"nb_epi"] <- p[i+1,1]
  tab[p[i,1],"poids_epis"] <- p[i+2,1]
}

# les épis ont été pesés avec leur sachet donc on enlève le poids du sachet
tab <- tab %>% mutate_all(.funs = as.numeric)
tab$poids_epis <- tab$poids_epis - 1475

# on met les poids neglligeables à 0 et on passe en gramme
tab$poids_epis <- ifelse(tab$poids_epis < 0 , 0 , tab$poids_epis/1000)


bac <- merge(bac , tab , by = "row.names" , all.x = T) %>% column_to_rownames(var = "Row.names")

save(bac , file = "bac")





# ajout variable bordure

rm(list = ls())

load("bac")

bac$bordure <- ifelse(bac$X == 1 | bac$X == 16 | bac$Y == 1 | bac$Y == 13 , "OUI" , "NON")

# verification
ggplot(bac , aes(x = X , y = Y , fill = bordure)) + geom_tile() + facet_wrap(~BAC)
# c'est ok

save(bac , file = "bac")

rm(list = ls())



# ajout covariables de position

load("bac")

traits <- c("N_flag","preco","hauteur","nb_epi","poids_epis")

for (t in traits){

  for (i in 1:nrow(bac)){
    
    # coordonnees de la plante focale
    xf <- bac[i,"X"]
    yf <- bac[i,"Y"]
    bacf <- bac[i,"BAC"]
    
    # precocite des plantes voisines
    p1 <- bac[which(bac$X == xf+1 & bac$Y == yf & bac$BAC == bacf),t]
    p2 <- bac[which(bac$X == xf-1 & bac$Y == yf & bac$BAC == bacf),t]
    p3 <- bac[which(bac$X == xf & bac$Y == yf+1 & bac$BAC == bacf),t]
    p4 <- bac[which(bac$X == xf & bac$Y == yf-1 & bac$BAC == bacf),t]
    p5 <- bac[which(bac$X == xf+1 & bac$Y == yf+1 & bac$BAC == bacf),t]
    p6 <- bac[which(bac$X == xf-1 & bac$Y == yf-1 & bac$BAC == bacf),t]
    p7 <- bac[which(bac$X == xf+1 & bac$Y == yf-1 & bac$BAC == bacf),t]
    p8 <- bac[which(bac$X == xf-1 & bac$Y == yf+1 & bac$BAC == bacf),t]
    
    # on garde que ceux qui sont pas vides
    la <- c()
    
    for (p in c(p1,p2,p3,p4,p5,p6,p7,p8)){
      if (is_empty(p)==F){la <- c(la,p)}
    }
    
    # On fait la moyenne et on la met dans le tableau
    bac[i,paste0(t,"_voisin")] <- mean(la , na.rm = T)
    
  }

}

rm(p1,p2,p3,p4,p5,p6,p7,p8,la,i,p,xf,yf,bacf,t,traits)

save(bac , file = "bac")





# Cr?ation d'une variable de presence/absence
bac$appel <- ifelse(is.na(bac$hauteur) == F , "present" , "absent")

save(bac , file = "bac")






# creation de variable nb_voisins

load("bac")

for (i in 1:nrow(bac)){
  
  # coordonnees de la plante focale
  xf <- bac[i,"X"]
  yf <- bac[i,"Y"]
  bacf <- bac[i,"BAC"]
  
  # precocite des plantes voisines
  p1 <- bac[which(bac$X == xf+1 & bac$Y == yf & bac$BAC == bacf),"appel"]
  p2 <- bac[which(bac$X == xf-1 & bac$Y == yf & bac$BAC == bacf),"appel"]
  p3 <- bac[which(bac$X == xf & bac$Y == yf+1 & bac$BAC == bacf),"appel"]
  p4 <- bac[which(bac$X == xf & bac$Y == yf-1 & bac$BAC == bacf),"appel"]
  p5 <- bac[which(bac$X == xf+1 & bac$Y == yf+1 & bac$BAC == bacf),"appel"]
  p6 <- bac[which(bac$X == xf-1 & bac$Y == yf-1 & bac$BAC == bacf),"appel"]
  p7 <- bac[which(bac$X == xf+1 & bac$Y == yf-1 & bac$BAC == bacf),"appel"]
  p8 <- bac[which(bac$X == xf-1 & bac$Y == yf+1 & bac$BAC == bacf),"appel"]
  
  nb <- 0
  
  for (p in c(p1,p2,p3,p4,p5,p6,p7,p8)){
    if (p == "present"){nb <- nb+1}
  }
  
  bac[i,"nb_voisin"] <- nb
  
  
}


save(bac , file = "bac")




# Ajout des donnees optomachine des grains recoltes

rm(list=ls())

# recuperation de tous les noms
file_names <- list.files(path = "./data_brute/opto_recolte_bac")

# vérif que tout s'est bien passe
# for (f in file_names){
#   tab <- read.table(paste0("./data_brute/opto_recolte_bac/",f,"/",f,"_METRO_CLASSIFICATION.dat") , header = T , sep = "\t" , dec = ",")
#   
#   a <- "S03" %in% tab$Réf..Ech
#   
#   if (a == T){print(f)}
# }
# ok y'a pas eu de cafouillage


# constitution d'un gros tableau avec toutes les donnees
library(data.table)

load("bac")

opto_recolte_bac <- data.frame()



for (f in file_names){
  
    tab <- read.table(paste0("./data_brute/opto_recolte_bac/",f,"/",f,"_METRO_CLASSIFICATION.dat") , header = T , sep = "\t" , dec = ",")
  
  # pour avoir le PMG
  tab2 <- fread(paste0("./data_brute/opto_recolte_bac/",f,"/",f,"_METRO_STATISTIQUES.dat") , header = T , sep = "\t" , dec = "," , nrows = 2)
  
  
  ind <- sapply(strsplit(f , "Mecarnot_") , "[" , 2)
  
  # cas du 201_1_2 qui est un rate pendant mesures
  if (ind == "201_1_2"){
    ind <- "201_1"
    tab$Réf..Ech <- "S02"
  }
  
  geno <- sapply(strsplit(ind , "_") , "[" , 1)
  grain <- sapply(strsplit(ind , "_") , "[" , 2)
  
  tab$ind <- ind
  tab$geno <- geno 
  tab$grain <- grain
  
  b <- bac[which(row.names(bac) == ind ) , c("BAC","bordure","semis","luz")]
  
  tab$BAC <- b$BAC
  tab$bordure <- b$bordure
  tab$semis <- b$semis
  tab$luz <- b$luz
  
  
  # sans classification cassés
  PMG <- 1000 * tab2$`Masse totale mesurée` / tab2$`Nb graines`
  tab$PMG <- PMG
  
  # Avec classification cassé
  tab$PMG2 <- tab2$`PMG Clas. C1`

  tab <- tab %>% relocate(luz , .before = Réf..Ech) %>% relocate(ind , .before = luz) %>% relocate(geno , .before = luz) %>% relocate(grain , .before = luz) %>% relocate(BAC , .before = luz) %>% relocate(bordure , .before = luz) %>% relocate(semis , .before = luz) %>% relocate(PMG , .before = luz) %>% relocate(PMG2 , .before = luz)
  
  opto_recolte_bac <- rbind(opto_recolte_bac , tab)
}


# renommage des variables comme pour opto

names(opto_recolte_bac)[10:108] <- c("Ref.Ech"  ,              
                 "Index",                     "Longueur"      ,       "Longueur.interieure" ,
                 "Largeur",              "Perimetre"      ,      "Perimetre.de.Crofton",
                 "Perimetre.convexe",    "Dimetre.Eq"     ,     "Surface"            ,
                 "Surface.convexe"  ,   "Surface.DiffCAP" ,    "Surface.DiffCEN"    ,
                 "Surface.DiffEAP"  ,   "Finesse"                ,   "Excentricite"             ,
                 "F.Feret"                ,  "Compacite"               ,  "Circularite"              ,
                 "Rugosite"                 , "Index.de.courbure",         "Moment.inertie.1"         ,
                 "Moment.inertie.2",          "Moment.inertie.3"  ,        "Symetrie"                 ,
                 "Moyenne.ndg"      ,         "Ecart.type.ndg"     ,       "Minimum.ndg"              ,
                 "Maximum.ndg"       ,        "Histo.Kurtose"      ,      "Histo.Moyenne"           ,
                 "Histo.Pic"         ,       "Histo.Assymetrie"    ,     "Histo.Ecart.type"        ,
                 "Histo.Variance"     ,      "Cooc.Uniformite"      ,    "Cooc.Contraste"          ,
                 "Cooc.Correlation"    ,     "Cooc.Variance.globale" ,   "Cooc.Homogeneïte"        ,
                 "Cooc.Somme.des.moyennes",  "Cooc.Somme.des.variances", "Cooc.Somme.des.entropies",
                 "Cooc.Entropie.globale",    "Cooc.Ecart.variance",      "Cooc.Ecart.entropie"     ,
                 "Cooc.Correlation.IC1"  ,   "Cooc.Correlation.IC2",     "RVB.R.Moy"                ,
                 "RVB.R.Ect"               ,  "RVB.R.Min"             ,    "RVB.R.Max"                ,
                 "RVB.V.Moy",                 "RVB.V.Ect"              ,   "RVB.V.Min"                ,
                 "RVB.V.Max" ,                "RVB.B.Moy"               ,  "RVB.B.Ect"                ,
                 "RVB.B.Min"  ,               "RVB.B.Max"                , "TSI.T.Moy"                ,
                 "TSI.T.Ect"   ,              "TSI.T.Min",                 "TSI.T.Max"                ,
                 "TSI.S.Moy"    ,             "TSI.S.Ect" ,                "TSI.S.Min"                ,
                 "TSI.S.Max"     ,            "TSI.I.Moy"  ,               "TSI.I.Ect"                ,
                 "TSI.I.Min"      ,           "TSI.I.Max"   ,              "CMJ.C.Moy"                ,
                 "CMJ.C.Ect"       ,          "CMJ.C.Min"    ,             "CMJ.C.Max"                ,
                 "CMJ.M.Moy"        ,         "CMJ.M.Ect"     ,            "CMJ.M.Min"                ,
                 "CMJ.M.Max"         ,        "CMJ.J.Moy"      ,           "CMJ.J.Ect"                ,
                 "CMJ.J.Min"          ,       "CMJ.J.Max"       ,          "Lab.L.Moy"                ,
                 "Lab.L.Ect"           ,      "Lab.L.Min"        ,         "Lab.L.Max"                ,
                 "Lab.a.Moy"            ,     "Lab.a.Ect"         ,        "Lab.a.Min"                ,
                 "Lab.a.Max"             ,    "Lab.b.Moy"          ,       "Lab.b.Ect"                ,
                 "Lab.b.Min"              ,   "Lab.b.Max" , "Classe" , "Masse.surfacique.gr.mm2" , "Poids.estime.g")




# sauvegarde de tout ça
save(opto_recolte_bac , file = "opto_recolte_bac")

rm(list = ls())




# ajout des variables au tableau bac
load("bac")
load("opto_recolte_bac")


# on ne garder que les epis du brin maitre et on prend des stats utiles

ajout_non_classe <- opto_recolte_bac %>% 
  mutate(count = 1) %>%
  filter(Ref.Ech == "S01") %>% 
  group_by(ind) %>% 
  summarise(surface_recolte_moy = mean(Surface) , 
            surface_recolte_min = min(Surface) , 
            surface_recolte_max = max(Surface) , 
            PMG = mean(PMG) , 
            poids_moy = mean(Poids.estime.g) , 
            poids_min = min(Poids.estime.g) , 
            poids_max = max(Poids.estime.g),
            nb_grain = sum(count)) %>% 
  column_to_rownames(var = "ind")

ajout_classe <- opto_recolte_bac %>% 
  mutate(count = 1) %>%
  filter(Ref.Ech == "S01" & Classe == 1) %>% 
  group_by(ind) %>% 
  summarise(surface_recolte_moy2 = mean(Surface) , 
            surface_recolte_min2 = min(Surface) , 
            surface_recolte_max2 = max(Surface) , 
            PMG2 = mean(PMG2) , 
            poids_moy2 = mean(Poids.estime.g) , 
            poids_min2 = min(Poids.estime.g) , 
            poids_max2 = max(Poids.estime.g)) %>%
  column_to_rownames(var = "ind")


bac <- merge(bac,ajout_non_classe , by = "row.names" , all.x = T) %>% 
  column_to_rownames(var = "Row.names") %>%
  merge(ajout_classe , by = "row.names" , all.x = T) %>%
  column_to_rownames(var = "Row.names")

save(bac , file = "bac")



# estimation des variances ------------------------------------------------
rm(list=ls())

load("opto_recolte")

# quels ind ont ete passes pour plusieurs epis ?
quel <- opto_recolte %>% filter(Ref.Ech == "S02")

quelind <- unique(quel$ind)

# on garde que ceux la
estim_var <- opto_recolte %>% filter(ind %in% quelind)




# Creation d'une variable epi digne de ce nom


estim_var[1,"epi"] <- epi <- 1

for (i in 2:nrow(estim_var)){
  
  if (estim_var[i-1,"Index"] > estim_var[i,"Index"]){epi <- epi+1}
  
  estim_var[i,"epi"] <- epi
}

estim_var <- estim_var %>% relocate(epi , .after = luz) %>% mutate_at(.vars = "epi" , .funs = as.factor)

save(estim_var , file = "estim_var")



rm()
# matrice genotypique et map ----------------------------------------------

rm(list = ls())

# on commence par matrice genotypique pour prediction genomique

load("Map_EPO_07_09_2021_Clement/SG_EPO_complet.Rdata")

load("bac")

g  <- unique(bac$geno)

row.names(SG) <- sapply(strsplit(row.names(SG) , split = "_") , "[" , 2)

genot <- SG[which(row.names(SG) %in% g),]

nrow(genot)
# il y a un probleme, il y a un genot en trop

sort(row.names(genot))
# le 27 apparait 2 fois

load("Map_EPO_07_09_2021_Clement/SG_EPO_complet.Rdata")
a <- genot[which(row.names(genot) == "27"),]
a <- rbind(a , SG[which(row.names(SG)=="EL4X_27"),])
# le 27.1 (le deuxieme 27 dans le tableau est le mauvais)

genot <- genot[-which(row.names(genot)=="27")[2],]


nrow(genot)

sort(row.names(genot))

a <- genot[which(row.names(genot) == "27"),]
a <- rbind(a , SG[which(row.names(SG)=="EL4X_27"),])
unique(a[1,]==a[2,])
# c'est good


class(genot[1,])
# mauvaise classe

genot_pred <- apply(genot , MARGIN = c(1,2) , FUN = as.numeric)

class(genot_pred[1,])
# bonne classe

# nettoyage

## deal with NAs
# nb NA
sum(is.na(genot_pred))

# fréquence des NA
sum(is.na(genot_pred)) / prod(dim(genot_pred))


# imputation des NA avec la classe génotypique la plus fréquente

genot2 <- apply(genot, 2, function(x){
  freq <- table(x)
  #  x[is.na(x)] <- as.integer(names(which.max(freq))) ne fonctionne pas ici car c'est codé en AB et pas 0 1 2
  x[is.na(x)] <- as.numeric(names(which.max(freq)))
  return(x)
})
sum(is.na(genot2))


class(genot2[1,])
#bon format

genot2 <- apply(genot2 , MARGIN = c(1,2) , FUN = as.numeric)

class(genot2[1,])
# bonne classe

dim(genot2)

# Filtre sur les maf
# frequences alléliques

p <- colMeans(genot2) / 2
q <- 1 - p


# MAF (minor allele frequency)

maf <- apply(cbind(p, q), 1, min)
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))


# Filtre sur MAF à 5%

sum(maf < 0.05)
genot.ok <- genot2[, maf >= 0.05]
dim(genot.ok)

class(genot.ok[1,])

# Vérification

p <- colMeans(genot.ok) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))



dim(genot.ok)

genot_pred <- genot.ok

save(genot_pred , file = "genot_pred")











# On passe aux donnees de map et matrice genotypique pour la gwas
rm(list = ls())

load("genot_pred")

map <- read.table("../donnees/Map_EPO_07_09_2021_Clement/Map_EPO_07_09_2021.csv" , sep = "," , dec = "." , header = T)

ncol(genot_pred)
nrow(map)
# Il y a plus de marqueurs dans genot que dans map. on ne garde que ceux de map


unique(snp %in% colnames(genot_pred))
# certains snp de map sont pas dans genot_pred

# on ne garde que le necessaire et on met tout dans le même ordre
snp <- which(map$Marker_ID %in% colnames(genot_pred))
row.names(map) <- map$Marker_ID

map <- map[snp,]

genot <- genot_pred[,row.names(map)]


ncol(genot)
unique(colnames(genot) == rownames(map))
# c'est good


map <- map %>% select(Marker_ID , Svevo_original_Start , Chr.Svevo_original) %>% rename(SNP = Marker_ID , pos = Svevo_original_Start , chr = Chr.Svevo_original)


unique(colnames(genot) == row.names(map))
# c'est good et dans le même ordre




dim(map)
dim(genot)


# On met map dans l'ordre
map <- map %>% arrange(chr,pos)

# On fait suivre pour genot
genot_gwas <- genot[,row.names(map)]


unique(colnames(genot_gwas) == row.names(map))
plot(match(row.names(map), colnames(genot_gwas)))
# c'est bon, tout est en ordre !



# ah je viens de voir ça : y'a des chromosomes "Un" et de NA :
which(map$chr == "Un")
which(is.na(map$chr) == T)

# on corrige ça :
map <- map %>% filter(chr != "Un" & is.na(chr)==F)

# dans map on fait une nouvel encodage des chr pour que ça passe dans certaines fct de manhatan plot

map$chr2 <- map$chr

for (i in 1:nrow(map)){
  if (map[i,"chr2"] == "1A"){map[i,"chr2"] <- 1}
  if (map[i,"chr2"] == "1B"){map[i,"chr2"] <- 2}
  if (map[i,"chr2"] == "2A"){map[i,"chr2"] <- 3}
  if (map[i,"chr2"] == "2B"){map[i,"chr2"] <- 4}
  if (map[i,"chr2"] == "3A"){map[i,"chr2"] <- 5}
  if (map[i,"chr2"] == "3B"){map[i,"chr2"] <- 6}
  if (map[i,"chr2"] == "4A"){map[i,"chr2"] <- 7}
  if (map[i,"chr2"] == "4B"){map[i,"chr2"] <- 8}
  if (map[i,"chr2"] == "5A"){map[i,"chr2"] <- 9}
  if (map[i,"chr2"] == "5B"){map[i,"chr2"] <- 10}
  if (map[i,"chr2"] == "6A"){map[i,"chr2"] <- 11}
  if (map[i,"chr2"] == "6B"){map[i,"chr2"] <- 12}
  if (map[i,"chr2"] == "7A"){map[i,"chr2"] <- 13}
  if (map[i,"chr2"] == "7B"){map[i,"chr2"] <- 14}
}



# On remet que les bons individus dans genot_gwas
genot_gwas <- genot_gwas[,row.names(map)]
dim(genot_gwas)
dim(map)

unique(colnames(genot_gwas) == row.names(map))
plot(match(row.names(map), colnames(genot_gwas)))
# c'est bon, tout est en ordre !

map$chr2 <- as.numeric(map$chr2)


save(genot_gwas , file = "genot_gwas")
save(map , file = "map")


rm(list=ls())







# donnees champ -----------------------------------------------------------

rm(list=ls())

tab <- read.table("data_brute/hauteur_champ.csv" , header = F , sep = ";" , dec = ".")


id <- tab[seq(1,nrow(tab),3),]

hauteur <- as.numeric(tab[seq(2,nrow(tab),3),])

rang <- as.numeric(tab[seq(3,nrow(tab),3),])


champ <- data.frame(ind = id , hauteur = hauteur , nb_epillets = 2*rang+1)


champ$passage <- sapply(strsplit(champ$ind , split = "_") , "[" , 1)

champ$planche <- sapply(strsplit(champ$ind , split = "_") , "[" , 2)

champ$selection <- sapply(strsplit(champ$ind , split = "_") , "[" , 3)


champ <- champ %>% column_to_rownames(var = "ind")

champ$parcelle <- paste0(champ$passage,champ$planche)

for (i in 1:nrow(champ)){
  if (champ[i,"parcelle"] == "11"){champ[i,"parcelle"] <- "1"}
  if (champ[i,"parcelle"] == "21"){champ[i,"parcelle"] <- "2"}
  if (champ[i,"parcelle"] == "31"){champ[i,"parcelle"] <- "3"}
  if (champ[i,"parcelle"] == "12"){champ[i,"parcelle"] <- "4"}
  if (champ[i,"parcelle"] == "22"){champ[i,"parcelle"] <- "5"}
  if (champ[i,"parcelle"] == "inconnu1inconnu1"){champ[i,"parcelle"] <- "6"}
  if (champ[i,"parcelle"] == "13"){champ[i,"parcelle"] <- "7"}
  if (champ[i,"parcelle"] == "inconnu2inconnu2"){champ[i,"parcelle"] <- "8"}
  if (champ[i,"parcelle"] == "33"){champ[i,"parcelle"] <- "9"}
  if (champ[i,"parcelle"] == "14"){champ[i,"parcelle"] <- "10"}
  if (champ[i,"parcelle"] == "24"){champ[i,"parcelle"] <- "11"}
  if (champ[i,"parcelle"] == "34"){champ[i,"parcelle"] <- "12"}
}


save(champ , file = "champ")




# formatage des données d'optomachine

rm(list = ls())

# correspondance entre les noms des etiquettes et le nom des id. l'info est contenue dans HAUTEUR_CHAMP_brute.xlsx

eti <- c("1Petit",
         "1gros",
         "1Non Trie",
         "2gros",
         "2Moyen",
         "3Non Trie",
         "1Moyen",
         "3Petit",
         "3Moyen",
         "3gros",
         "2Non Trie",
         "2Petit")

id <- c("1_1_Petit",
        "3_3_gros",
        "inconnu2_inconnu2_Non Trie",
        "3_1_gros",
        "2_1_Moyen",
        "inconnu1_inconnu1_Non Trie",
        "1_2_Moyen",
        "1_3_Petit",
        "3_4_Moyen",
        "2_2_gros",
        "1_4_Non Trie",
        "2_4_Petit")

corre <- data.frame(etiquette = eti , id = id)


file_names <- list.files(path = "./data_brute/opto_recolte_champ")

opto_recolte_champ <- data.frame()

for (f in file_names){
  
  tab <- read.table(paste0("./data_brute/opto_recolte_champ/",f,"/",f,"_METRO_CLASSIFICATION.dat") , header = T , sep = "\t" , dec = ",")
  
  # pour avoir le PMG
  tab2 <- fread(paste0("./data_brute/opto_recolte_champ/",f,"/",f,"_METRO_STATISTIQUES.dat") , header = T , sep = "\t" , dec = "," , nrows = 2)
  
  
  ind <- sapply(strsplit(f , "Rep ") , "[" , 2)
  
  eti <- sapply(strsplit(ind , "_") , "[" , 1)
  pp <- corre[which(corre$etiquette == eti),"id"]
  epi <- sapply(strsplit(ind , "_") , "[" , 2)
  
  id <- paste0(pp,"_",epi)
  
  pas <- sapply(strsplit(id , "_") , "[" , 1)
  
  pl <- sapply(strsplit(id , "_") , "[" , 2)
  
  tab$id <- id
  tab$passage <- pas
  tab$planche <- pl
  
  # sans classification cassés
  PMG <- 1000 * tab2$`Masse totale mesurée` / tab2$`Nb graines`
  tab$PMG <- PMG
  
  # Avec classification cassé
  tab$PMG2 <- tab2$`PMG Clas. C1`
  
  opto_recolte_champ <- rbind(opto_recolte_champ , tab)
}

names(opto_recolte_champ)[1:99] <- c("Ref.Ech"  ,              
                 "Index",                     "Longueur"      ,       "Longueur.interieure" ,
                 "Largeur",              "Perimetre"      ,      "Perimetre.de.Crofton",
                 "Perimetre.convexe",    "Dimetre.Eq"     ,     "Surface"            ,
                 "Surface.convexe"  ,   "Surface.DiffCAP" ,    "Surface.DiffCEN"    ,
                 "Surface.DiffEAP"  ,   "Finesse"                ,   "Excentricite"             ,
                 "F.Feret"                ,  "Compacite"               ,  "Circularite"              ,
                 "Rugosite"                 , "Index.de.courbure",         "Moment.inertie.1"         ,
                 "Moment.inertie.2",          "Moment.inertie.3"  ,        "Symetrie"                 ,
                 "Moyenne.ndg"      ,         "Ecart.type.ndg"     ,       "Minimum.ndg"              ,
                 "Maximum.ndg"       ,        "Histo.Kurtose"      ,      "Histo.Moyenne"           ,
                 "Histo.Pic"         ,       "Histo.Assymetrie"    ,     "Histo.Ecart.type"        ,
                 "Histo.Variance"     ,      "Cooc.Uniformite"      ,    "Cooc.Contraste"          ,
                 "Cooc.Correlation"    ,     "Cooc.Variance.globale" ,   "Cooc.Homogeneïte"        ,
                 "Cooc.Somme.des.moyennes",  "Cooc.Somme.des.variances", "Cooc.Somme.des.entropies",
                 "Cooc.Entropie.globale",    "Cooc.Ecart.variance",      "Cooc.Ecart.entropie"     ,
                 "Cooc.Correlation.IC1"  ,   "Cooc.Correlation.IC2",     "RVB.R.Moy"                ,
                 "RVB.R.Ect"               ,  "RVB.R.Min"             ,    "RVB.R.Max"                ,
                 "RVB.V.Moy",                 "RVB.V.Ect"              ,   "RVB.V.Min"                ,
                 "RVB.V.Max" ,                "RVB.B.Moy"               ,  "RVB.B.Ect"                ,
                 "RVB.B.Min"  ,               "RVB.B.Max"                , "TSI.T.Moy"                ,
                 "TSI.T.Ect"   ,              "TSI.T.Min",                 "TSI.T.Max"                ,
                 "TSI.S.Moy"    ,             "TSI.S.Ect" ,                "TSI.S.Min"                ,
                 "TSI.S.Max"     ,            "TSI.I.Moy"  ,               "TSI.I.Ect"                ,
                 "TSI.I.Min"      ,           "TSI.I.Max"   ,              "CMJ.C.Moy"                ,
                 "CMJ.C.Ect"       ,          "CMJ.C.Min"    ,             "CMJ.C.Max"                ,
                 "CMJ.M.Moy"        ,         "CMJ.M.Ect"     ,            "CMJ.M.Min"                ,
                 "CMJ.M.Max"         ,        "CMJ.J.Moy"      ,           "CMJ.J.Ect"                ,
                 "CMJ.J.Min"          ,       "CMJ.J.Max"       ,          "Lab.L.Moy"                ,
                 "Lab.L.Ect"           ,      "Lab.L.Min"        ,         "Lab.L.Max"                ,
                 "Lab.a.Moy"            ,     "Lab.a.Ect"         ,        "Lab.a.Min"                ,
                 "Lab.a.Max"             ,    "Lab.b.Moy"          ,       "Lab.b.Ect"                ,
                 "Lab.b.Min"              ,   "Lab.b.Max" , "Classe" , "Masse.surfacique.gr.mm2" , "Poids.estime.g")


# creation d'une variable parcelle


opto_recolte_champ$parcelle <- paste0(opto_recolte_champ$passage,opto_recolte_champ$planche)

for (i in 1:nrow(opto_recolte_champ)){
  if (opto_recolte_champ[i,"parcelle"] == "11"){opto_recolte_champ[i,"parcelle"] <- "1"}
  if (opto_recolte_champ[i,"parcelle"] == "21"){opto_recolte_champ[i,"parcelle"] <- "2"}
  if (opto_recolte_champ[i,"parcelle"] == "31"){opto_recolte_champ[i,"parcelle"] <- "3"}
  if (opto_recolte_champ[i,"parcelle"] == "12"){opto_recolte_champ[i,"parcelle"] <- "4"}
  if (opto_recolte_champ[i,"parcelle"] == "22"){opto_recolte_champ[i,"parcelle"] <- "5"}
  if (opto_recolte_champ[i,"parcelle"] == "inconnu1inconnu1"){opto_recolte_champ[i,"parcelle"] <- "6"}
  if (opto_recolte_champ[i,"parcelle"] == "13"){opto_recolte_champ[i,"parcelle"] <- "7"}
  if (opto_recolte_champ[i,"parcelle"] == "inconnu2inconnu2"){opto_recolte_champ[i,"parcelle"] <- "8"}
  if (opto_recolte_champ[i,"parcelle"] == "33"){opto_recolte_champ[i,"parcelle"] <- "9"}
  if (opto_recolte_champ[i,"parcelle"] == "14"){opto_recolte_champ[i,"parcelle"] <- "10"}
  if (opto_recolte_champ[i,"parcelle"] == "24"){opto_recolte_champ[i,"parcelle"] <- "11"}
  if (opto_recolte_champ[i,"parcelle"] == "34"){opto_recolte_champ[i,"parcelle"] <- "12"}
}


opto_recolte_champ <- opto_recolte_champ %>% 
  relocate(PMG2 , .before = Ref.Ech) %>% 
  relocate(id , .before = PMG2) %>%
  relocate(parcelle , .before = PMG2) %>%
  relocate(passage , .before = PMG2) %>%
  relocate(planche , .before = PMG2) %>%
  relocate(PMG , .before = PMG2)


save(opto_recolte_champ , file = "opto_recolte_champ")






# ajout des variables au tableau champ

rm(list=ls())
load("champ")
load("opto_recolte_champ")



ajout_non_classe <- opto_recolte_champ %>% 
  mutate(count = 1) %>%
  group_by(id) %>% 
  summarise(surface_recolte_moy = mean(Surface) , 
            surface_recolte_min = min(Surface) , 
            surface_recolte_max = max(Surface) , 
            PMG = mean(PMG) , 
            poids_moy = mean(Poids.estime.g) , 
            poids_min = min(Poids.estime.g) , 
            poids_max = max(Poids.estime.g),
            nb_grain = sum(count)) %>% 
  column_to_rownames(var = "id")

ajout_classe <- opto_recolte_champ %>% 
  mutate(count = 1) %>%
  filter(Classe == 1) %>% 
  group_by(id) %>% 
  summarise(surface_recolte_moy2 = mean(Surface) , 
            surface_recolte_min2 = min(Surface) , 
            surface_recolte_max2 = max(Surface) , 
            PMG2 = mean(PMG2) , 
            poids_moy2 = mean(Poids.estime.g) , 
            poids_min2 = min(Poids.estime.g) , 
            poids_max2 = max(Poids.estime.g)) %>%
  column_to_rownames(var = "id")



champ <- merge(champ,ajout_non_classe , by = "row.names" , all.x = T) %>% 
  column_to_rownames(var = "Row.names") %>%
  merge(ajout_classe , by = "row.names" , all.x = T) %>%
  column_to_rownames(var = "Row.names")



save(champ , file = "champ")

