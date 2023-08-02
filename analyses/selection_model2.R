rm(list=ls())

setwd("~/Stage/Analyses")

library(ggplot2)
library(tidyverse)


load("../donnees/bac")



# graphiques --------------------------------------------------------------

don <- bac %>% filter(geno !="INCONNU" & appel == "present") #%>% mutate(nb_voisin = as.factor(nb_voisin))

ggplot(don , aes(x=PMG , y=prot_recolte)) + geom_point()
ggplot(don , aes(x=PMG , y=hauteur)) + geom_point()
ggplot(don , aes(x=PMG , y=nb_epi)) + geom_point()
ggplot(don , aes(x=PMG , y=poids_epis)) + geom_point()

ggplot(don , aes(y=hauteur , x = hauteur_voisin)) + geom_point()
ggplot(don , aes(y=PMG , x = hauteur_voisin)) + geom_point()
ggplot(don , aes(y=nb_epi , x = hauteur_voisin)) + geom_point()
ggplot(don , aes(y=poids_epis , x = hauteur_voisin)) + geom_point()
ggplot(don , aes(y=surface_recolte_moy , x = hauteur_voisin)) + geom_point()

ggplot(don , aes(y=surface_recolte_moy , x=hauteur)) + geom_point()
ggplot(don , aes(y=surface_recolte_moy , x=hauteur_voisin)) + geom_point()
ggplot(don , aes(y=surface_recolte_moy , x=nb_epi)) + geom_point()
ggplot(don , aes(y=surface_recolte_moy , x=poids_epis)) + geom_point()


ggplot(don , aes(y=hauteur , x = nb_voisin)) + geom_boxplot()
ggplot(don , aes(y=PMG , x = nb_voisin)) + geom_boxplot()
ggplot(don , aes(y=nb_epi , x = nb_voisin)) + geom_boxplot()
ggplot(don , aes(y=poids_epis , x = nb_voisin)) + geom_boxplot()
ggplot(don , aes(y=surface_recolte_moy , x = nb_voisin)) + geom_boxplot()

ggplot(don , aes(x=hauteur_voisin , y=poids_epis_voisin)) + geom_point()
ggplot(don , aes(x=hauteur_voisin , y=hauteur)) + geom_point()
ggplot(don , aes(x=hauteur , y=poids_epis_voisin)) + geom_point()



# model complet

v <- "PMG"
f <- as.formula(paste(v,"~ BAC + semis + bordure + N_flag_voisin + hauteur_voisin + preco_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin"))
step(lm(f , data=don), direction="both")


v <- "surface_recolte_moy"
f <- as.formula(paste(v,"~ BAC + semis + bordure + N_flag_voisin + hauteur_voisin + preco_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin"))
step(lm(f , data=don), direction="both")


v <- "hauteur"
f <- as.formula(paste(v,"~ BAC + semis + bordure + N_flag_voisin + hauteur_voisin + preco_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin"))
step(lm(f , data=don), direction="both")



v <- "poids_epis"
f <- as.formula(paste(v,"~ BAC + semis + N_flag_voisin + hauteur_voisin + preco_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin"))
step(lm(f , data=don), direction="both")



don <- bac %>% filter(geno != "INCONNU" & appel == "present" & hauteur > 50 & nb_grain > 10)

v <- "prot_recolte"
f <- as.formula(paste(v,"~ BAC + semis + bordure + N_flag_voisin + hauteur_voisin + preco_voisin + nb_epi_voisin + poids_epis_voisin + nb_voisin"))
step(lm(f , data=don), direction="both")
