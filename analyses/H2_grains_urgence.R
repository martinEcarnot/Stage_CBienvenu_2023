rm(list = ls())




setwd("~/Stage/donnees")
#setwd("~/Documents/INRA/SelPhen_2023/Stage_CBienvenu_2023/donnees/")

library(tidyverse)
library(lme4)
library(lmerTest)


load("opto_recolte_bac")


don <- opto_recolte_bac %>% filter(Classe == 1 & semis == "06/01" & geno != "INCONNU")

don$count <- 1

tmp <- don %>% group_by(ind) %>% summarise(NGE = sum(count))


mod <- lmer(Surface ~ (1|geno) + BAC, data = don)
Vg <- VarCorr(mod)$geno[1]
Vr <- (sigma(mod))^2
Vg/(Vg + Vr)



don_epi <- don %>% mutate(BAC = as.numeric(BAC) , geno = as.numeric(geno)) %>% group_by(ind) %>% summarise(Surface = mean(Surface) , geno = mean(geno) , BAC = mean(BAC)) %>% mutate(BAC = as.factor(BAC),geno = as.factor(geno))



mod <- lmer(Surface ~ (1|geno) + BAC, data = don_epi)
Vg <- VarCorr(mod)$geno[1]
#Vb <- VarCorr(mod)$BAC[1]
Vr <- (sigma(mod))^2
Vg/(Vg + Vr)










load("opto")

mod <- lmer(prot_semis ~ (1|geno) , data = opto)
Vg <- VarCorr(mod)$geno[1]
Vr <- (sigma(mod))^2
Vg/(Vg + Vr)



mod <- lmer(Surface ~ (1|geno) , data = opto)
Vg <- VarCorr(mod)$geno[1]
Vr <- (sigma(mod))^2
Vg/(Vg + Vr)
