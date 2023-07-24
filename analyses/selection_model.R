rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)
library(ggpubr)


load("../donnees/bac")



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




# l'ordi qui boss ---------------------------------------------------------

don <- bac %>% filter(geno != "INCONNU" & appel == "present")


# hauteur
formule1 <- hauteur ~ geno + semis + luz + BAC + bordure + luz/BAC + BAC/bordure + hauteur_voisin
formule2 <- hauteur ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + hauteur_voisin
formule3 <- hauteur ~ geno + semis + BAC + bordure + BAC:bordure + hauteur_voisin
formule4 <- hauteur ~ geno + semis + luz + bordure + luz:bordure + hauteur_voisin


mod <- lm(formule1 <- hauteur ~ geno + semis + luz + BAC + bordure + luz/BAC + BAC/bordure + hauteur_voisin , data = don)
drop1(mod , .~. , test = "F")

mod <- lm(formule1 <- hauteur ~ BAC, data = don)
drop1(mod , .~. , test = "F")

lmer(hauteur ~ (1|geno) + luz + hauteur_voisin , data = don)

step(lm(formule1 , data=don), direction="both")
# hauteur ~ geno + semis + bordure + hauteur_voisin    AIC : 5082.13

step(lm(formule2 , data=don), direction="both")
# hauteur ~ geno + semis + hauteur_voisin + BAC:bordure    AIC : 5094.39

step(lm(formule3 , data=don), direction="both")
# hauteur ~ geno + semis + bordure + hauteur_voisin  AIC=5082.13

step(lm(formule4 , data=don), direction="both")
# hauteur ~ geno + semis + bordure + hauteur_voisin   AIC=5082.13


# preco
formule1 <- preco ~ geno + semis + BAC + luz + bordure + luz:BAC + BAC:bordure + preco_voisin
formule2 <- preco ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + preco_voisin
formule3 <- preco ~ geno + semis + BAC + bordure + BAC:bordure + preco_voisin
formule4 <- preco ~ geno + semis + luz + bordure + luz:bordure + preco_voisin


step(lm(formule1 , data=don), direction="both")
# preco ~ geno + semis + BAC + bordure + preco_voisin     AIC = 7374.77

step(lm(formule2 , data=don), direction="both")
# preco ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + preco_voisin   AIC : 7382.73

step(lm(formule3 , data=don), direction="both")
# preco ~ geno + semis + BAC + bordure + preco_voisin  AIC=7374.77

step(lm(formule4 , data=don), direction="both")
# preco ~ geno + semis + luz + bordure + preco_voisin  AIC=7408.31


# N_flag
formule1 <- N_flag ~ geno + semis + BAC + luz + bordure + luz:BAC + BAC:bordure + N_flag_voisin
formule2 <- N_flag ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + N_flag_voisin
formule3 <- N_flag ~ geno + semis + BAC + bordure + BAC:bordure + N_flag_voisin
formule4 <- N_flag ~ geno + semis + luz + bordure + luz:bordure + N_flag_voisin

step(lm(formule1 , data=don), direction="both")
# N_flag ~ semis + BAC + bordure + N_flag_voisin  AIC=-2038.71

step(lm(formule2 , data=don), direction="both")
# N_flag ~ semis + luz + N_flag_voisin + luz:BAC + BAC:bordure  AIC=-2034.97

step(lm(formule3 , data=don), direction="both")
# N_flag ~ semis + BAC + bordure + N_flag_voisin  AIC=-2038.71

step(lm(formule4 , data=don), direction="both")
# N_flag ~ semis + luz + bordure + N_flag_voisin  AIC=-2001.32


# nb_epi
formule1 <- nb_epi ~ geno + semis + BAC + luz + bordure + luz:BAC + BAC:bordure + nb_epi_voisin
formule2 <- nb_epi ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + nb_epi_voisin
formule3 <- nb_epi ~ geno + semis + BAC + bordure + BAC:bordure + nb_epi_voisin
formule4 <- nb_epi ~ geno + semis + luz + bordure + luz:bordure + nb_epi_voisin

step(lm(formule1 , data=don), direction="both")
# nb_epi ~ semis + BAC + bordure + nb_epi_voisin   AIC=-135.58

step(lm(formule2 , data=don), direction="both")
#nb_epi ~ semis + luz + nb_epi_voisin + luz:BAC + BAC:bordure   AIC=-130.83 

step(lm(formule3 , data=don), direction="both")
# nb_epi ~ semis + BAC + bordure + nb_epi_voisin   AIC=-135.58

step(lm(formule4 , data=don), direction="both")
# nb_epi ~ semis + luz + bordure + nb_epi_voisin   AIC=-122.55


# poids_epis
formule1 <- poids_epis ~ geno + semis + BAC + luz + bordure + luz:BAC + BAC:bordure + poids_epis_voisin + nb_epi
formule2 <- poids_epis ~ geno + semis + luz + BAC %in% luz + bordure %in% BAC + poids_epis_voisin + nb_epi
formule3 <- poids_epis ~ geno + semis + BAC + bordure + BAC:bordure + poids_epis_voisin + nb_epi
formule4 <- poids_epis ~ geno + semis + luz + bordure + luz:bordure + poids_epis_voisin + nb_epi

step(lm(formule1 , data=don), direction="both")
# poids_epis ~ semis + BAC + nb_epi  AIC=1194.26

step(lm(formule2 , data=don), direction="both")
# poids_epis ~ semis + luz + nb_epi + luz:BAC   AIC=1194.26

step(lm(formule3 , data=don), direction="both")
# poids_epis ~ semis + BAC + nb_epi   AIC=1194.26

step(lm(formule4 , data=don), direction="both")
#  poids_epis ~ semis + luz + nb_epi   AIC=1198.19




# On conserve les models realistes avec les plus petits AIC et on regarde si c'est valide et les significativites

# hauteur
mod <- lm(hauteur ~ geno + semis + bordure + hauteur_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")

# preco
mod <- lm(preco ~ geno + semis + BAC + bordure + preco_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")

# N_flag
mod <- lm(N_flag ~ semis + BAC + bordure + N_flag_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")

# nb_epi
mod <- lm(nb_epi ~ semis + BAC + bordure + nb_epi_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")
# pas bon

# poids_epis
mod <- lm(poids_epis ~ semis + BAC + nb_epi , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")
# pas bon





# tests pour nb voisin

mod <- lm(hauteur ~ geno + semis + bordure + hauteur_voisin + nb_voisin , data = don)
summary(mod)
drop1(mod , .~. , test = "F")



mod <- lm(preco ~ geno + semis + BAC + bordure + preco_voisin + nb_voisin , data = don)
summary(mod)
drop1(mod , .~. , test = "F")


mod <- lm(N_flag ~ semis + BAC + bordure + N_flag_voisin + nb_voisin , data = don)
summary(mod)
drop1(mod , .~. , test = "F")
# SIGNIFICATIF



# nb_epi
mod <- lm(nb_epi ~ semis + BAC + bordure + nb_epi_voisin + nb_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")
# pas bon

# poids_epis
mod <- lm(poids_epis ~ semis + BAC + nb_epi + nb_voisin , data = don)
par(mfrow = c(2,2))
plot(mod)
drop1(mod , .~. , test = "F")
# pas bon












# effet du génotype sur les traits des voisin -----------------------------

don <- bac %>% filter(geno != "INCONNU" & appel == "present")

# hauteur
formule <- hauteur_voisin ~ geno + semis + bordure

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared


# preco
formule <- preco_voisin ~ geno + semis + BAC + bordure

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared


# N_flag
formule <- N_flag_voisin ~ geno + semis + BAC + bordure

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared


# nb_epi
formule <- nb_epi_voisin ~ geno + semis + BAC + bordure

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared



# poids_epis
formule <- poids_epis_voisin ~ geno + semis + BAC + nb_epi

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared



# nb_voisin
formule <- nb_voisin ~ geno + semis + bordure

mod <- lm(formule , data = don)

drop1(mod , .~. , test = "F")
summary(mod)$r.squared
