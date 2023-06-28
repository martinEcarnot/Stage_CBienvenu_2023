setwd("~/APIMET/UE_2_Biodiversite_innovation_varietale_et_production_de_semences/2.3_Genetique_et_creation_varietale/Cartographie et QTL/Rapport évalué/Partie Jacques")

rm(list = ls())

library(anyLib)
anyLib(c("apercu", "corpcor", "data.table", "mlmm", "MM4LMM", "qqman" , "tidyverse" , "ggplot2" , "rrBLUP"))

# Importation des données -------------------------------------------------

genot <- read.table("../marqueurs_DIC2_SILUR.csv" , sep = ";" , dec = "." , header = T , na.strings = "-") %>% rename(Ind_id = geno)
genot.mat <- as.matrix(genot[, 2:ncol(genot)])
rownames(genot.mat) <- genot$Ind_id


y <- read.table("../pheno_champ.csv" , sep = ";" , dec = "." , header = T) %>% filter(is.na(epiaison)==F) %>% group_by(geno) %>% summarise(epiaison = mean(epiaison)) %>% column_to_rownames(var = "geno")


map <- read.table("../map_DIC2_SILUR.csv" , sep = ";" , dec = "." , header = T) %>% rename(Chr = chromo_map , Pos = position_map , SNP = marker)

ggplot(y , aes(x = epiaison)) + geom_histogram(bins = 20)




# Deal with NAs -----------------------------------------------------------

# nb NA
sum(is.na(genot.mat))

# fréquence des NA
sum(is.na(genot.mat)) / prod(dim(genot.mat))


# imputation des NA avec la classe génotypique la plus fréquente

genot.i <- apply(genot.mat, 2, function(x){
  freq <- table(x)
#  x[is.na(x)] <- as.integer(names(which.max(freq))) ne fonctionne pas ici car c'est codé en AB et pas 0 1 2
  x[is.na(x)] <- names(which.max(freq))
  return(x)
})
#sum(is.na(genot.i))





# Recodage de A et B en 0 et 2 --------------------------------------------

genot.imp <- apply(genot.i, c(1,2), function(x){
  if(x == "A"){x <- 0}
  if(x == "B"){x <- 2}
  return(x)
})



# On retire les maf < 5% --------------------------------------------------

# frequences alléliques

p <- colMeans(genot.imp) / 2
q <- 1 - p


# MAF (minor allele frequency)

maf <- apply(cbind(p, q), 1, min)
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))


# Filtre sur MAF à 5%

sum(maf < 0.05)
genot.ok <- genot.imp[, maf >= 0.05]
dim(genot.ok)


# Vérification

p <- colMeans(genot.ok) / 2
q <- 1 - p
maf <- apply(cbind(p, q), 1, min)
hist(maf, col = "grey", main = "", breaks = 50, xlim = c(0, 0.5))


# Nettoyage de l'environnement

rm(genot, genot.mat, genot.imp, maf, p, q , genot.i)





# Carte -------------------------------------------------------------------

# Filtre sur les SNP conservés :
map <- map[map$SNP %in% colnames(genot.ok), ]


# Tri par chromosome et position :
map <- map[order(map$Pos), ]
map <- map[order(map$Chr), ]

# Tri de y pour que ce soit dans le même ordre que map
y2 <- y
y <- data.frame()

for (i in row.names(genot.ok)){
  y[i,"epiaison"] <- y2[i,"epiaison"]
}
rm(y2)

## Vérification du tri des données
# *Note* : pour s'assurer que toutes les données sont bien dans le même ordre

plot(match(rownames(genot.ok), rownames(y)))
plot(match(map$SNP, colnames(genot.ok)))


# c'est ok




# matrice de kinship ------------------------------------------------------



# Matrice d'apparentement génomique :

## standardisation (centage-réduction) de la matrice de génotypage

p <- colMeans(genot.ok) / 2
q <- 1 - p
genot.scaled <- scale(genot.ok, center = 2 * p, scale = sqrt(2 * p * q))

# calcul de la matrice d'apparentement

# *Note* : `tcrossprod(X)` est quivalent à `X %*% t(X)`

K <- tcrossprod(genot.scaled) / ncol(genot.scaled)


# vérification que la matrice est définie positive

is.positive.definite(K)
K <- make.positive.definite(K)
is.positive.definite(K)


# representation graphique
image(K)





# GWAS --------------------------------------------------------------------

## Modèle simple locus avec MM4LMM

### Inférence


mmest <- MMEst(Y = y$epiaison, X = genot.ok, VarList = list(Additive = K, Error = diag(length(y$epiaison))))
out.test <- AnovaTest(mmest)
res.mm4lmm <- cbind(map,P = sapply(out.test, function(x){x["Xeffect", "pval"]}))


### Représentation graphique

#### QQ-plot

qq(res.mm4lmm$P)



# Transformation des chromosomes 1A 1B etc... en 1 2 3 psq sinon la fonction manhatan plot marche pas

res.mm4lmm$Chr2 <- res.mm4lmm$Chr

#### Manhattan plot

#### Manhattan plot
ggplot(res.mm4lmm , aes(x = Pos , y = -log10(P))) + geom_point() + facet_grid(.~Chr2) + geom_hline(yintercept = -log10(0.05/nrow(res.mm4lmm)) , col = "red") + theme(axis.text.x = element_blank() , strip.switch.pad.grid = unit(0,"cm"))


for (i in 1:nrow(res.mm4lmm)){
  if (res.mm4lmm[i,"Chr2"] == "1A"){res.mm4lmm[i,"Chr"] <- 1}
  if (res.mm4lmm[i,"Chr2"] == "1B"){res.mm4lmm[i,"Chr"] <- 2}
  if (res.mm4lmm[i,"Chr2"] == "2A"){res.mm4lmm[i,"Chr"] <- 3}
  if (res.mm4lmm[i,"Chr2"] == "2B"){res.mm4lmm[i,"Chr"] <- 4}
  if (res.mm4lmm[i,"Chr2"] == "3A"){res.mm4lmm[i,"Chr"] <- 5}
  if (res.mm4lmm[i,"Chr2"] == "3B"){res.mm4lmm[i,"Chr"] <- 6}
  if (res.mm4lmm[i,"Chr2"] == "4A"){res.mm4lmm[i,"Chr"] <- 7}
  if (res.mm4lmm[i,"Chr2"] == "4B"){res.mm4lmm[i,"Chr"] <- 8}
  if (res.mm4lmm[i,"Chr2"] == "5A"){res.mm4lmm[i,"Chr"] <- 9}
  if (res.mm4lmm[i,"Chr2"] == "5B"){res.mm4lmm[i,"Chr"] <- 10}
  if (res.mm4lmm[i,"Chr2"] == "6A"){res.mm4lmm[i,"Chr"] <- 11}
  if (res.mm4lmm[i,"Chr2"] == "6B"){res.mm4lmm[i,"Chr"] <- 12}
  if (res.mm4lmm[i,"Chr2"] == "7A"){res.mm4lmm[i,"Chr"] <- 13}
  if (res.mm4lmm[i,"Chr2"] == "7B"){res.mm4lmm[i,"Chr"] <- 14}
}

res.mm4lmm$Chr <- as.numeric(res.mm4lmm$Chr)

manhattan(res.mm4lmm, chr = "Chr", bp = "Pos", p = "P",
          suggestiveline = FALSE,
          genomewideline = -log10(0.05/nrow(res.mm4lmm)))


### Effets des SNP significatifs

signifSNPs <- subset(res.mm4lmm, P <= 0.05/nrow(res.mm4lmm))
sapply(mmest[signifSNPs$SNP], function(x) x$Beta)[2, ]

### Illustration des effets

for (i in signifSNPs$SNP) {
  boxplot(y$epiaison ~ genot.ok[, i], varwidth = T, main = i, xlab = "")
  mtext(text = c("n = ", table(genot.ok[, i])), side = 1,
        at = c(0.5, 1, 2), line = 3)
}

## Modèle multi-locus avec MLMM

### Inférence

# *Notes* :
#   
#   - on limite le nombre de steps à 6
# 
# - l'argument `nbchunks` permet de limiter l'utilisation de la RAM

map$Chr2 <- map$Chr
for (i in 1:nrow(map)){
  if (map[i,"Chr2"] == "1A"){map[i,"Chr"] <- 1}
  if (map[i,"Chr2"] == "1B"){map[i,"Chr"] <- 2}
  if (map[i,"Chr2"] == "2A"){map[i,"Chr"] <- 3}
  if (map[i,"Chr2"] == "2B"){map[i,"Chr"] <- 4}
  if (map[i,"Chr2"] == "3A"){map[i,"Chr"] <- 5}
  if (map[i,"Chr2"] == "3B"){map[i,"Chr"] <- 6}
  if (map[i,"Chr2"] == "4A"){map[i,"Chr"] <- 7}
  if (map[i,"Chr2"] == "4B"){map[i,"Chr"] <- 8}
  if (map[i,"Chr2"] == "5A"){map[i,"Chr"] <- 9}
  if (map[i,"Chr2"] == "5B"){map[i,"Chr"] <- 10}
  if (map[i,"Chr2"] == "6A"){map[i,"Chr"] <- 11}
  if (map[i,"Chr2"] == "6B"){map[i,"Chr"] <- 12}
  if (map[i,"Chr2"] == "7A"){map[i,"Chr"] <- 13}
  if (map[i,"Chr2"] == "7B"){map[i,"Chr"] <- 14}
}

map$Chr <- as.numeric(map$Chr)

mygwas <- mlmm(Y = y$epiaison , X = genot.ok, K = K, maxsteps = 6, nbchunks = 2)

mygwas$step_table

### Manhattan plots

#### Step 1 : pas de cofacteur

plot_fwd_GWAS(x = mygwas, step = 1, snp_info = map, pval_filt = 0.1)

#### Step 2 : 1 cofacteur

plot_fwd_GWAS(x = mygwas, step = 2, snp_info = map, pval_filt = 0.1)

#### Step 3 : 2 cofacteurs

plot_fwd_GWAS(x = mygwas, step = 3, snp_info = map, pval_filt = 0.1)

#### Meilleur modèle selon le critère "mbonf"

plot_step_table(mygwas, "maxpval")

plot_opt_GWAS(x = mygwas, opt = "mbonf", snp_info = map, pval_filt = 0.1)

### QQplots

#### Step 1

qqplot_fwd_GWAS(x = mygwas, nsteps = 6)

#### Meilleur modèle selon le critère "mbonf"

qqplot_opt_GWAS(x = mygwas, opt = "mbonf")

### Partition de variance

plot_step_RSS(mygwas)


### Coefficients et p-valeurs des cofacteurs

mygwas$opt_mbonf$coef

### Illustration des effets

for (i in mygwas$opt_mbonf$cof) {
  boxplot(y$epiaison ~ genot.ok[, i], varwidth = T, main = i, xlab = "")
  mtext(text = c("n = ", table(genot.ok[, i])), side = 1,
        at = c(0.5, 1, 2), line = 3)
}







# Prédiction génomique ----------------------------------------------------
library(rrBLUP)

rm(list = ls())


genot <- as.matrix(read.table("../marqueurs_DIC2_SILUR.csv" , sep = ";" , dec = "." , header = T , na.strings = "-") %>% column_to_rownames(var = "geno"))

y2 <- as.matrix(read.table("../pheno_champ.csv" , sep = ";" , dec = "." , header = T) %>% filter(is.na(epiaison)==F & env != "mauguio2010" & env != "cappelle2013") %>% group_by(geno) %>% summarise(epiaison = mean(epiaison)) %>% column_to_rownames(var = "geno"))


# On met tout dans le même ordre

Pheno <- data.frame()

for (i in row.names(genot)){
  Pheno[i,"epiaison"] <- y2[i,"epiaison"]
}
rm(y2)


plot(match(rownames(genot), rownames(Pheno)))
# c'est ok




# Conversion des A et B en -1 et 1

Markers <- apply(genot, c(1,2), function(x){
  if(is.na(x) == T){return(x)}
  if(x == "A"){x <- -1}
  if(x == "B"){x <- 1}
  return(x)
})

rm(genot)

Markers <- apply(Markers , 2 , as.numeric)

Pheno <- as.matrix(Pheno)

#remove markers with more than 50% missing data
col <- c()
for (i in ncol(Markers)){
  if( sum(is.na(Markers[,i])) / nrow(Markers) >= 0.5){col <- c(col , i)}
}

#Markers <- Markers[,-col] # là pas besoin car col = Null donc pas de plus de 50% de données manquantes
rm(col)

#what if markers are NA?
#impute with A.mat
impute <- A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute2 <- impute$imputed

rm(impute , Markers , i)



# Cross validation :

trait <- "epiaison"

#prop_train <- 0.7

cycles <- 100

accuracy  <-  data.frame(prop = 0 , accuracy = 0)

res <- data.frame(predit = 0 , obs = 0 , prop = 0)

for (p in seq(0.2 , 0.8 , 0.05)){
  prop_train <- p

  for(r in 1:cycles){
    train <- sample(1:nrow(Pheno), floor(prop_train*nrow(Pheno)))
    test <- setdiff(1:nrow(Pheno),train)
    Pheno_train <- Pheno[train,]
    m_train <- Markers_impute2[train,]
    Pheno_valid <- Pheno[test,]
    m_valid <- Markers_impute2[test,]
    
  #  yield=(Pheno_train[,trait])
    
    answer <- mixed.solve(Pheno_train , Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
    u <- as.matrix(answer$u)
    pred_valid <-  m_valid %*% u
    
    accuracy <- rbind(accuracy , c(p , cor(pred_valid, Pheno_valid , use="complete")) )
    
    pred <- pred_valid + rep(answer$beta , length(pred_valid))
    
    res <- rbind(res , data.frame(predit = pred , obs = Pheno_valid , prop = rep(p , length(Pheno_valid))))
  }
}

res <- res[-1,]
accuracy <- accuracy[-1,]
accuracy$prop <- as.factor(accuracy$prop)
res$prop <- as.factor(res$prop)

rm(answer , m_train , m_valid , pred_valid , train , test , r , Pheno_train , u , pred , Pheno_valid)


ggplot(as.data.frame(accuracy) , aes(y = accuracy , x = prop)) + geom_boxplot() + labs(title = "Boxplot des valeurs d'accuracy sur 100 cycles de cross validation" , y = "Accuracy") 

ggplot(subset(res , accuracy == "0.7")  , aes(x = predit , y = obs)) + geom_point() + labs(title = "Valeurs prédites VS observées pour toutes les validations" , x = "Valeurs prédites" , y = "Valeurs observées")

#ggplot(res , aes(x = predit , y = obs)) + geom_point() + labs(title = "Valeurs prédites VS observées pour toutes les validations" , x = "Valeurs prédites" , y = "Valeurs observées") + facet_wrap(~prop)

cor(subset(res , prop == 0.7)$predit , subset(res , prop == 0.7)$obs)**2

rm(res,accuracy , cycles , p , trait , prop_train)







# Calcul de progres genetique ---------------------------------------------

rm(list = ls())


# on enlève mauguio2010 et cappelle2013 pour meilleure prédiction
genot <- as.matrix(read.table("../marqueurs_DIC2_SILUR.csv" , sep = ";" , dec = "." , header = T , na.strings = "-") %>% column_to_rownames(var = "geno"))

y2 <- as.matrix(read.table("../pheno_champ.csv" , sep = ";" , dec = "." , header = T) %>% filter(is.na(epiaison)==F & env != "mauguio2010" & env != "cappelle2013") %>% group_by(geno) %>% summarise(epiaison = mean(epiaison)) %>% column_to_rownames(var = "geno"))


# On met tout dans le même ordre

Pheno <- data.frame()

for (i in row.names(genot)){
  Pheno[i,"epiaison"] <- y2[i,"epiaison"]
}
rm(y2)


plot(match(rownames(genot), rownames(Pheno)))
# c'est ok




# Conversion des A et B en -1 et 1

Markers <- apply(genot, c(1,2), function(x){
  if(is.na(x) == T){return(x)}
  if(x == "A"){x <- -1}
  if(x == "B"){x <- 1}
  return(x)
})

rm(genot)

Markers <- apply(Markers , 2 , as.numeric)



#what if markers are NA?
#impute with A.mat
impute <- A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers <- impute$imputed

rm(impute , i)


row.names(Markers) <- row.names(Pheno)

Markers <- as.data.frame(Markers)



### Calcul de S
# Sélection des individus étant 1 pour les SNPs 1572 1573 et 1574 (1 = B = 2 pour la GWAS = précocifiant)

mark_sel <- subset(Markers , SNP1572 == 1 & SNP1573 == 1 & SNP1574 == 1)

phen_sel <- as.data.frame(Pheno[which(row.names(Pheno) %in% row.names(mark_sel)),])
names(phen_sel) <- "epiaison"
row.names(phen_sel) <- row.names(mark_sel)


# rpz graphique
S <- mean(phen_sel$epiaison)

ggplot(Pheno , aes(x = epiaison)) + geom_histogram(bins = 15) + geom_vline(xintercept = mean(Pheno$epiaison) , col = "blue") + geom_vline(xintercept = S , col = "red") + annotate("text" , label = "new µ" , x = S-2 , y = 1 , col = "red" , size = 5) + annotate("text" , label = "µ" , x = mean(Pheno$epiaison) + 1 , y = 1 , col = "blue" , size = 5) + labs(y = "Count" , title = "Histogramme des valeurs d'épiaison")

ggplot(phen_sel , aes(x = epiaison)) + geom_histogram(bins = 15) + labs(y = "Count" , title = "Histogramme des valeurs d'épiaison des individus sélectionnés")

### Calcul du h2

# On crée des enfants : on fait les 2 parmi n croisements possibles
# pour chaque enfant, la valeur du marqueur est la moyenne des marqueurs des 2 parents (les deux parents sont homozygotes)

# enfants <- Markers[1,]
# 
# moy_p <- c()
# 
# for (i in 70:nrow(Markers)){
#   for (j in c((i+1) : nrow(Markers))){
#     if (i < nrow(Markers) - 1){
#       parents <- Markers[c(i,j),]
#       enfants <- rbind(enfants , apply(parents , 2 , mean) )
#       moy_p <- c(moy_p , (Pheno[i,"epiaison"] + Pheno[j,"epiaison"])/2)
#     }
#   }
#   print(i)
# }
# 
# enfants <- enfants[-1,]
# 
# save(enfants , file = "enfants")
# save(moy_p , file = "moy_p")


#rm(enfants,enf,moy_p)

load("enfants")
load("moy_p")


# établissement du modèle :

answer <- mixed.solve(as.matrix(Pheno) , Z=as.matrix(Markers), K=NULL, SE = FALSE, return.Hinv=FALSE)
u <- as.matrix(answer$u)


# Prédiction du phénotype des enfants :

pred_enf <- as.matrix(enfants) %*% u + rep(answer$beta , nrow(enfants))

library(ggpubr)
ggplot(as.data.frame(cbind(pred_enf,moy_p)) , aes(x = moy_p , y = V1)) + geom_point() + geom_smooth(method = "lm" , col = "red" , se = F) + stat_regline_equation(size = 5) + labs(title = "Phenotype prédit des enfants VS moyenne de leurs parents" , y = "Phénotypes prédits des enfants" , x = "Moyenne phénotypique des parents")

plot(moy_p , pred_enf)


h2 <- lm(pred_enf ~ moy_p)$coefficients[2]

prog <- h2*(S - mean(Pheno$epiaison))
