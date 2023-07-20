rm(list=ls())

setwd("~/Stage/Analyses")

load("../donnees/champ")

library(tidyverse)




don <- champ %>% filter(!passage %in% c("inconnu1","inconnu2") & !planche %in% c("inconnu1","inconnu2"))


moy <- don %>% group_by(passage,planche) %>% summarise_at(.vars = c("hauteur" , "nb_epillets") , .funs = mean , na.rm = T)


ggplot(moy , aes(x = passage , y = planche , fill = hauteur)) + geom_tile()

ggplot(moy , aes(x = passage , y = planche , fill = nb_epillets)) + geom_tile()


ggplot(champ , aes(x = hauteur , y = nb_epillets)) + geom_point() + geom_smooth(method = "lm" , col = "red")


mod <- lm(hauteur ~ selection , data = champ)
summary(mod)
drop1(mod , .~. , test = "F")


mod <- lm(nb_epillets ~ selection , data = champ)
summary(mod)
drop1(mod , .~. , test = "F")
