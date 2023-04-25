rm(list = ls())

# setwd("~/Stage/donnees")
setwd("~/Documents/INRA/SelPhen_2023/Stage_CBienvenu_2023/donnees/")

library(nirsextra)
library(tidyverse)
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





# Donnees optomachine -----------------------------------------------------

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
  # On commence par creer une variable de presence ou non dans les bacs pour indiquer si ça a germe. Initialisation a OUI et toutes les fois ou c'est pas le cas en vrai, on modifie en NON :
  
  bac$germi <- "OUI"
  
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
  
  
  # on modifie bac : on met "NON" dans la colonne germi pour dire que ça a pas germe et on ajoute une ligne avec l'individu qui le remplace 
  c <- 1 # compteur utile dans le if des genotypes inconnus
  
  for (i in 1:nrow(correc)){
    ind <- correc[i,"ind"]
    new_ind <- correc[i,"new_ind"]
    
    # Si on sait quel genotype a remplace : on duplique la ligne de base et on change son numero de grain, puis on change la valeur de germination de la ligne de base
    if (is.na(correc[ind , "grain"]) == FALSE){
      bac[new_ind,] <- bac[ind,]
      bac[new_ind,"grain"] <- correc[ind,"grain"]
      bac[ind , "germi"] <- "NON"
    }
    
    # Si on ne sait pas quel genotype a ete seme (les NA dans correc), on change juste le nom du genotype en INCONNU
    
    if (is.na(correc[ind , "grain"]) == TRUE){
      bac[ind,"geno"] <- "INCONNU"
      bac[ind,"grain"] <- c
      c <- c+1
    }
  }

  rm(c)
# On actualise les rownames
row.names(bac) <- paste0(bac$geno , "_" , bac$grain)
  
  rm(ind , new_ind , i)
  # ça c'etait pour toutes les donnees surlignees en orange sur le papier sans annotations roses (sauf pour la boite 12 où c'est pris en compte malgre les annotations roses).
  
  # il y a d'autres donnees plus chiantes qu'on peut pas simplement remplacer comme ça car les boites et les rangs entre l'individu remplace et le remplaçant ne correspondent pas. Du coup faut aller chercher à la main à quels genotype ils correspondent. Ce sont toutes les donnees entourrees en rose sur le format papier (sauf celles de la boite 12)
  
  {
    # recuperation des infos sur le genotype remplaçant
    a <- subset(bac , boite == 6 & Rang.OK == "E" & Debut == 3 & grain == 3)
    
    # genotype remplace
    ind <- "305_2"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  
  {
    a <- subset(bac , boite == 10 & Rang.OK == "F" & Debut == 3 & grain == 3)
    ind <- "488_7"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  {
    a <- subset(bac , boite == 10 & Rang.OK == "C" & Debut == 10 & grain == 4)
    ind <- "188_2"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  {
    a <- subset(bac , boite == 2 & Rang.OK == "C" & Debut == 8 & grain == 2)
    ind <- "65_2"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  
  
  {
    a <- subset(bac , boite == 9 & Rang.OK == "D" & Debut == 7 & grain == 1)
    ind <- "428_5"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  {
    a <- subset(bac , boite == 8 & Rang.OK == "C" & Debut == 2 & grain == 2)
    ind <- "428_6"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  
  {
    a <- subset(bac , boite == 7 & Rang.OK == "H" & Debut == 2 & grain == 2)
    ind <- "311_6"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  
  {
    a <- subset(bac , boite == 4 & Rang.OK == "A" & Debut == 1 & grain == 1)
    ind <- "489_6"
    
    new_geno <- a$geno
    new_grain <- a$grain + 6
    new_ind <- paste0(new_geno,"_",new_grain)
    
    bac[new_ind , ] <- bac[ind , ]
    bac[ind,"germi"] <- "NON"
    
    bac[new_ind,"geno"] <- new_geno
    bac[new_ind,"grain"] <- new_grain
  }
  
  
  rm(ind,new_geno,new_ind,new_grain,a,correc)
  
  # On enleve les variables inutiles 
  
  bac$X.1 <- bac$Fin <- bac$boite <- bac$boite.1 <- bac$Debut <- bac$RANG <- bac$Rang.OK <- bac$graines <- bac$gen <- NULL
  
  save(bac , file = "bac")
  
  rm(bac)
  
}

# Ajout des données de germination du 20/03/2023

# load("bac")

# On trouve les genotypes qui ont pas germe le 20/03 : les y et les x correspondent a ce qu'il y a sur papier. On precise germi = "OUI" pour bien cibler les genotypes qui ont remplace ceux qui ont pas germe a la premiere fois (les genotypes remplancants ont ete note "OUI" par defaut)

# no <- c(which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 1 & bac$X == 10),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 1 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 2 & bac$X == 5),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 2 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 4 & bac$X == 11),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 5 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 7 & bac$X == 12),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 8 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 10 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 13 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 1 & bac$Y == 13 & bac$X == 14),
#         
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 3 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 4 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 5 & bac$X == 1),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 6 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 8 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 8 & bac$X == 13),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 9 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 2 & bac$Y == 13 & bac$X == 9),
#         
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 1 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 2 & bac$X == 5),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 5 & bac$X == 2),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 5 & bac$X == 3),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 5 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 6 & bac$X == 1),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 6 & bac$X == 10),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 6 & bac$X == 11),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 8 & bac$X == 5),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 8 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 9 & bac$X == 6),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 9 & bac$X == 10),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 10 & bac$X == 7),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 10 & bac$X == 8),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 11 & bac$X == 11),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 11 & bac$X == 12),
#         which(bac$germi == "OUI" & bac$BAC == 3 & bac$Y == 11 & bac$X == 15),
#         
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 1 & bac$X == 15),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 2 & bac$X == 12),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 3 & bac$X == 3),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 7 & bac$X == 1),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 7 & bac$X == 9),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 11 & bac$X == 5),
#         which(bac$germi == "OUI" & bac$BAC == 4 & bac$Y == 12 & bac$X == 16),
#         
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 1 & bac$X == 1),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 4 & bac$X == 8),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 7 & bac$X == 4),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 7 & bac$X == 6),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 8 & bac$X == 1),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 10 & bac$X == 16),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 11 & bac$X == 12),
#         which(bac$germi == "OUI" & bac$BAC == 5 & bac$Y == 13 & bac$X == 2),
#         
#         which(bac$germi == "OUI" & bac$BAC == 6 & bac$Y == 3 & bac$X == 15),
#         which(bac$germi == "OUI" & bac$BAC == 6 & bac$Y == 6 & bac$X == 8),
#         which(bac$germi == "OUI" & bac$BAC == 6 & bac$Y == 13 & bac$X == 8)
# )
# 
# bac[no,"germi"] <- "NON"
# 
# save(bac , file = "bac")
# 
# 
# rm(list = ls())
# }
# fonction pour passer une matrice comme une colonne dans un dataframe : sp2df (il faut que les psectres soient sous format matrice)