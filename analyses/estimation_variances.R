rm(list=ls())

setwd("~/Stage/Analyses")

library(tidyverse)
library(lme4)
library(lmerTest)

load("../donnees/estim_var")

names(estim_var)
unique(estim_var$semis)

#############
# LE MODEL ##
#############

don <- estim_var %>% filter(semis == "06/01")

mod <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

mod






mod <- lmer(Surface ~ (1|geno) + (1|BAC) + (1|geno:(ind:epi)) , data = don)

mod





mod <- lmer(Surface ~ semis + (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = estim_var)
anova(mod)

mod <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = estim_var)
mod

mod <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = estim_var %>% filter(semis == "06/01"))
mod






don <- estim_var %>% filter(semis == "06/01")

mod1 <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

mod2 <- lmer(Surface ~ (1|geno) + (1|geno:ind) , data = don)

anova(mod1,mod2)





mod1 <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

mod2 <- lmer(Surface ~ (1|geno) + (1|geno:(ind:epi)) , data = don)

anova(mod1,mod2)




mod1 <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

mod2 <- lmer(Surface ~ (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

anova(mod1,mod2)


# tout semble significatif









mod1 <- lmer(Surface ~ (1|geno) + (1|geno:ind) + (1|geno:(ind:epi)) , data = don)

mod2 <- lmer(Surface ~ (1|geno) , data = don)

mod1
mod2


anova(mod1,mod2)
