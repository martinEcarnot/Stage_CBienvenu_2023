rm(list=ls())

setwd("~/APIMET/UE_1_Outils_de_diagnostic_pour_l_innovation_varietale/ECUE_1.3_Methodes_et_outils_de_modelisation_des_cultures")

library(tidyverse)
library(ggplot2)

sensitivity_eb = read.csv("resultats/sensitivity_eb.csv", header = T, sep=";") %>% mutate(sensi = "eb" , WUE = Rdt / eau_conso)

sensitivity_HI = read.csv("resultats/sensitivity_HI.csv", header = T, sep =";") %>% mutate(sensi = "HI" , WUE = Rdt / eau_conso)

sensitivity_interception= read.csv("resultats/sensitivity_interception.csv", header = T, sep=";") %>% mutate(sensi = "interception" , WUE = Rdt / eau_conso)

sensitivity_phenologie = read.csv("resultats/sensitivity_phenologie.csv", header = T, sep=";") %>% mutate(sensi = "phenologie" , WUE = Rdt / eau_conso)

sensitivity_Prof_Rac = read.csv("resultats/sensitivity_Prof_Rac.csv", header = T, sep=";") %>% mutate(sensi = "Prof_Rac" , WUE = Rdt / eau_conso)

sensitivity_resp_WS = read.csv("resultats/sensitivity_resp_WS.csv", header = T, sep=";") %>% mutate(sensi = "resp_WS" , WUE = Rdt / eau_conso)

sensitivity_SF = read.csv("resultats/sensitivity_SF.csv", header = T, sep=";") %>% mutate(sensi = "SF" , WUE = Rdt / eau_conso)


grp <- read.table("Groupe_annee.csv" , sep = ";" , header = T)

sensitivity <- rbind(sensitivity_eb , sensitivity_HI , sensitivity_interception , sensitivity_phenologie , sensitivity_Prof_Rac , sensitivity_resp_WS , sensitivity_SF) %>% merge(grp , by = "Annee")

rm(sensitivity_eb , sensitivity_HI , sensitivity_interception , sensitivity_phenologie , sensitivity_Prof_Rac , sensitivity_resp_WS , sensitivity_SF , grp)


#### 1. script à completer en intégrant Prof_Rac, resp_WS et SF
####  2. Refaire script pour sensibilité mutli-critère (Rdt_potentiel, conso_eau, WUE)
#### 3. Analyser le cv pour les différents types de climats




CV <- sensitivity %>% group_by(Annee , sensi , groupe_annee) %>% summarise(CV_rdt = sd(Rdt)/mean(Rdt) , CV_rdt_pot = sd(Rdt_pot)/mean(Rdt_pot) , CV_wue = sd(WUE)/mean(WUE)) %>% filter(Annee != "1992") %>% mutate_at(.vars = "groupe_annee" , .funs = as.factor)

ggplot(CV , aes(x = Annee , y = CV_rdt)) + geom_col() + facet_grid(groupe_annee~sensi)
ggplot(CV , aes(x = sensi , y = CV_rdt , fill = groupe_annee)) + geom_boxplot(alpha = 0.5)

ggplot(CV , aes(x = Annee , y = CV_rdt_pot)) + geom_col() + facet_grid(~sensi)
ggplot(CV , aes(x = sensi , y = CV_rdt_pot)) + geom_boxplot()

ggplot(CV , aes(x = Annee , y = CV_wue)) + geom_col() + facet_grid(~sensi)
ggplot(CV , aes(x = sensi , y = CV_wue)) + geom_boxplot()








'
# SCRIPT DE BASE :

par(mfrow = c(2,2))

sens_eb = aggregate(sensitivity_eb$Rdt, by = list(sensitivity_eb$Annee), sd) / aggregate(sensitivity_eb$Rdt, by = list(sensitivity_eb$Annee), mean)
sens_eb[,1] = unique(sensitivity_eb$Annee)
colnames(sens_eb) = c("Annee","Cv")
barplot(sens_eb$Cv, names.arg= unique(sensitivity_eb$Annee), ylim=c(0,0.1), main ="sensitivity_eb", xlab = "Annee", ylab = "Cv")

sens_phenologie = aggregate(sensitivity_phenologie$Rdt, by = list(sensitivity_phenologie$Annee), sd) / aggregate(sensitivity_phenologie$Rdt, by = list(sensitivity_phenologie$Annee), mean)
sens_phenologie[,1] = unique(sensitivity_eb$Annee)
colnames(sens_phenologie) = c("Annee","Cv")
barplot(sens_phenologie$Cv, names.arg= unique(sens_phenologie$Annee), ylim=c(0,0.1), main = "sensitivity_phenologie", xlab = "Annee", ylab = "Cv")


sens_interception = aggregate(sensitivity_interception$Rdt, by = list(sensitivity_interception$Annee), sd) / aggregate(sensitivity_interception$Rdt, by = list(sensitivity_interception$Annee), mean)
sens_interception[,1] = unique(sensitivity_interception$Annee)
colnames(sens_interception) = c("Annee","Cv")
barplot(sens_interception$Cv, names.arg= unique(sens_interception$Annee), ylim=c(0,0.1), main = "sensitivity_interception", xlab = "Annee", ylab = "Cv")


sens_eb$var = c("eb")
sens_phenologie$var = c("pheno")
sens_interception$var = c("interception")

final_data= rbind(sens_eb, sens_phenologie, sens_interception)

par(mfrow = c(1,1))
boxplot(final_data$Cv ~ final_data$var, xlab = "variable_type", ylab = "Cv")
'


