rm(list = ls())
setwd("~/Stage/analyses")

library(anyLib)
anyLib(c("apercu", "corpcor", "data.table", "mlmm", "MM4LMM", "qqman" , "tidyverse" , "ggplot2" , "rrBLUP" , "lme4"))

# Importation des données -------------------------------------------------

load("../donnees/bac")

load("../donnees/map")

load("../donnees/genot_gwas")

bac <- bac %>% filter(geno != "INCONNU")



# parametres pour extraction des BLUPs ------------------------------------

# donnees utilisees
don <- bac %>% filter(is.na(hauteur)==F)

# modele
f <- preco ~ (1|geno) + BAC + semis


{

# GWAS --------------------------------------------------------------------

# Extraction des BLUPs
mod <- lmer(f, data = don)
y <- ranef(mod)$geno %>% rename(blup = "(Intercept)")
rm(mod,don)

# Matrice de kinship
genot.ok <- genot_gwas[row.names(y),]
p <- colMeans(genot.ok) / 2
q <- 1 - p
genot.scaled <- scale(genot.ok, center = 2 * p, scale = sqrt(2 * p * q))
K <- tcrossprod(genot.scaled) / ncol(genot.scaled)
K <- make.positive.definite(K)
rm(genot.scaled,p,q)


## Modèle simple locus avec MM4LMM

### Inférence

mmest <- MMEst(Y = y$blup, X = genot.ok, VarList = list(Additive = K, Error = diag(length(y$blup))))
out.test <- AnovaTest(mmest)
res.mm4lmm <- cbind(map,P = sapply(out.test, function(x){x["Xeffect", "pval"]}))


### Représentation graphique

#### QQ-plot

qq(res.mm4lmm$P)


#### Manhattan plot
manhattan(res.mm4lmm, chr = "chr2", bp = "pos", p = "P",
          suggestiveline = FALSE,
          genomewideline = -log10(0.05/nrow(res.mm4lmm)))


### Effets des SNP significatifs

signifSNPs <- subset(res.mm4lmm, P <= 0.05/nrow(res.mm4lmm))
sapply(mmest[signifSNPs$SNP], function(x) x$Beta)[2, ]

### Illustration des effets

for (i in signifSNPs$SNP) {
  boxplot(y$blup ~ genot.ok[, i], varwidth = T, main = i, xlab = "")
  mtext(text = c("n = ", table(genot.ok[, i])), side = 1,
        at = c(0.5, 1, 2), line = 3)
}


}
