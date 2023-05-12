rm(list=ls())

setwd("~/Stage/Analyses")



library(prospectr)
library(signal)
library(rrBLUP)
library(tidyverse)

rm(list=ls())

load("../donnees/spectres")
load("../donnees/Map_EPO_07_09_2021_Clement/SG_EPO_complet.Rdata")

# format of data ----------------------------------------------------------

row.names(SG) <- sapply(strsplit(row.names(SG) , split = "_") , "[",2)

lignes <- row.names(spectres$brutes) %>% strsplit(split = "_") %>% sapply("[",1) %>% unique() %>% sort()

lignes <- lignes[! lignes %in% setdiff(lignes,row.names(SG))]

geno <- SG[lignes,] %>% as.data.frame() %>% apply(MARGIN = 1 , FUN = as.numeric) #la fonction apply avec margin = 1 sort une transposée

# verification que geno et SG contiennent la même info au même endroit 
# a <- t(geno) == SG[lignes,]
# unique(a[1,])
# ok

row.names(geno) <- colnames(SG)

colnames(geno) <- lignes

rm(SG)

geno <- geno %>% na.omit()

# Tout est dans le bon ordre normalement

# Filters on genomic data and computation of kinship matrix ---------------

# Filter on MAF
###############################

geno <- geno/2

p <- rowMeans(geno) # average frequency of the reference allele
summary(p)

ToRemove <- which(p <= 0.025 | p >= 0.975) # Remove markers with a MAF below 5%
length(ToRemove)

geno <- geno[-ToRemove, ]
dim(geno)

# Compute Kinship (matA1)
###############################

p <- rowMeans(geno)
q <- 1-p

genot.ok <- 2*t(geno)
rm(geno)

genot.scaled <- scale(genot.ok, center=2*p, scale=sqrt(4*p*q))

matA1 <- tcrossprod(genot.scaled) / ncol(genot.scaled)    # tcrossprod(X) is equivalent to X %*% t(X)
# matA1 is your genomic kinship

rm(p, q, genot.scaled, ToRemove)




# Are the spectra under genetic determinism ? -----------------------------

# Fit a GBLUP model to each wavelength to estimate the genomic and residual variances
# along the spectrum
###############################

spec <- spectres$brutes 

spec$geno <- g <- sapply(strsplit(row.names(spec) , split = "_"),"[",1)

spec <- spec[sort(row.names(spec)),] %>% filter(geno %in% lignes) %>% mutate(geno = NULL)

# verification que tout est dans le bon ordre
unique(row.names(matA1) == sort(row.names(matA1)))
unique(row.names(spec) == sort(row.names(spec)))

# matrice de design
design <- matrix(ncol = nrow(matA1) , nrow = nrow(spec) , data = 0)


j <- 1

for (i in 1:nrow(design)){
  design[i,j] <- 1
  
  if (i < nrow(design) & g[i] != g[i+1]){
    j <- j+1
  }
}



# variances de chaque longueur d'onde (TRES TRES TRES LONG)
GenomicVariance <- ResidualVariance <- rep(NA, nrow(spec))

for (i in 1:nrow(spec)) {
  print(i)
  
  mod4 <- mixed.solve(y = spec[,i], K = matA1 , Z = design)
  GenomicVariance[i] <- mod4$Vu
  ResidualVariance[i] <- mod4$Ve
  rm(mod4)
}

PropGenomicVariance <- GenomicVariance/(GenomicVariance+ResidualVariance)*100

#Graphical representation of the proportion of variance explained by genomics along the spectrum

par(mar=c(4, 4, 4, 4))
plot(seq(400, 2498, by=2), PropGenomicVariance, type="l", xlab="lambda (nm)",
     xlim=c(400, 3000), ylim=c(0, 100), ylab="Proportion of of variance explained by genomics")
polygon(c(400, seq(400, 2498, by=2), 2498), c(0, PropGenomicVariance, 0), col = "brown1")
polygon(c(400, seq(400, 2498, by=2), 2498), c(100, PropGenomicVariance, 100), col = "dodgerblue4")
legend(2500, 90, c("Residual", "Genomic"), lty = 0, bty = "n", fill = c("dodgerblue4", "brown1"), cex=1)

rm(g,GenomicVariance,i,j,PropGenomicVariance,ResidualVariance,design,spec)



# phenomique a proprement parler ------------------------------------------

rm(list = ls())
load("../donnees/spectres")
load("../donnees/opto")



# Here we use GBLUP but of course you could use any GS model

spct <- spectres$brutes         # Choose a pretreatment
spct2 <- scale(spct, center=T, scale=T) # scale absorbance at each wavelength (predictor)
matH <- tcrossprod(spct2)/ncol(spct2)     # Compute the hyperspectral similarity matrix

Nenvt=8   # Number of environments
Nind=nrow(pheno) # Number of varieties

Nrep=25   # Number of repetition for the cross validation
Nout=30   # Number of varieties in the predicted set


AccuHBLUP <- AccuGBLUP <- matrix(NA, Nrep, Nenvt)
colnames(AccuHBLUP) <- colnames(AccuGBLUP) <- colnames(pheno)[2:ncol(pheno)]


for (envt in 2:9) {
  print(envt)
  phenotype <- pheno[, envt]
  
  for (rep in 1:Nrep) {
    
    valid <- sample(Nind, Nout)
    phenoTrain <- phenotype
    phenoTrain[valid] <- NA
    
    gblup <- mixed.solve(y=phenoTrain, K=matA1)
    hblup <- mixed.solve(y=phenoTrain, K=matH)
    AccuGBLUP[rep,(envt-1)] <- cor(gblup$u[valid], phenotype[valid], use="complete.obs")
    AccuHBLUP[rep,(envt-1)] <- cor(hblup$u[valid], phenotype[valid], use="complete.obs")
  }
}

# Graphical representation of the results
###############################

# Boxplots

par(mfrow=c(1, 2), mar=c(8, 4, 2, 2))
# Predictive abilities obtained with GBLUP (genomic predictions)
boxplot(AccuGBLUP, ylab="Predictive abilities",
        ylim=c(-0.2, 1), las=2, col=c("blue", rep("lightblue", 7)),
        main="GBLUP (genomic prediction)", cex.axis=1)
abline(v=1.5)

# Predictive abilities obtained with GBLUP (genomic predictions)
boxplot(AccuHBLUP, ylab="Predictive abilities",
        ylim=c(-0.2, 1), las=2, col=c("red", rep("indianred", 7)),
        main="HBLUP (phenomic prediction)", cex.axis=1)
abline(v=1.5)

# Scatter plot HBLUP vs GBLUP

par(mfrow=c(1, 1), mar=c(4, 4, 4, 4))
plot(colMeans(AccuGBLUP), colMeans(AccuHBLUP),
     xlim=c(0, 1), ylim=c(0, 1),
     xlab="Predictive ability GBLUP (Genomic selection)",
     ylab="Predictive ability HBLUP (Phenomic selection)")
abline(a=0,b=1)
points(colMeans(AccuGBLUP)[1], colMeans(AccuHBLUP)[1],
       pch=22, col="red", bg="red") # highlight the reference environment
legend("bottomright", col=c("red", "black"),
       legend=c("Reference environment (with NIRS)","Other environments (without NIRS)"),
       pch=c(22, 1), cex=1)


# end of script
par(mfrow=c(1,2))
image(matA1 , main = "matA1")
image(matH , main = "matH")
