rm(list=ls())

setwd("~/Stage/Analyses")

set.seed(123)

load("../donnees/bac")

library(tidyverse)
#library(multcomp)

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


load("../donnees/opto")

pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F & appel == "present" & semis == "06/01")


moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

moy_geno <- moy_geno[which(moy_geno$geno %in% unique(pop$geno)),]

rm(opto)




# traits dont on regarde l'évolution

# ggplot(pop , aes(sample = preco)) + geom_qq() + geom_qq_line(col = "red")
# ggplot(popsel , aes(sample = preco)) + geom_qq() + geom_qq_line(col = "red")
# 
# ggplot(pop , aes(sample = hauteur)) + geom_qq() + geom_qq_line(col = "red")
# ggplot(popsel , aes(sample = hauteur)) + geom_qq() + geom_qq_line(col = "red")
# 
# ggplot(pop , aes(sample = N_flag)) + geom_qq() + geom_qq_line(col = "red")
# ggplot(popsel , aes(sample = N_flag)) + geom_qq() + geom_qq_line(col = "red")
# 
# ggplot(pop , aes(sample = nb_epi)) + geom_qq() + geom_qq_line(col = "red")
# ggplot(popsel , aes(sample = nb_epi)) + geom_qq() + geom_qq_line(col = "red")
# 
# ggplot(pop , aes(sample = poids_epis)) + geom_qq() + geom_qq_line(col = "red")
# ggplot(popsel , aes(sample = poids_epis)) + geom_qq() + geom_qq_line(col = "red")




traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis")


# Sélection light pour bonne estimation -----------------------------------

# sur grain individuel

don <- data.frame()
test <- data.frame()
nb_geno <- data.frame()
i <- 1

for (p in c(0.5,0.45,0.4,0.35,0.3,0.25,0.2)){
  
  a <- quantile(pop$Surface , probs = 1-p)
  ind <- pop %>% filter(Surface > a)
  
  gen <- moy_geno[which(moy_geno$Surface > quantile(moy_geno$Surface , probs = 1-p)) , "geno"]
  lot <- pop %>% filter(geno %in% gen$geno)
  rm(gen,a)
  
  alea <- pop[sample(nrow(pop) , round(p*nrow(pop))),]
  
  par_bac <- pop %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  
  tmp <- ind %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_ind")] <- tmp[,traits]
  
  tmp <- lot %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_lot")] <- tmp[,traits]
  
  par_bac <- as.data.frame(par_bac)

  par_bac$p <- p
  par_bac$n_ind <- nrow(ind)
  par_bac$n_lot <- nrow(lot)
  par_bac$i <- i_p(p)

  don <- rbind(don,par_bac)
  
  tmp <- data.frame(ind = length(unique(ind$geno)) , lot = length(unique(lot$geno)) , alea = length(unique(alea$geno)) , p = p , i = i_p(p))
  
  nb_geno <- rbind(nb_geno , tmp)
  
  for (t in c(traits)){
    
    alt <- "less"
    if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,t])){alt <- "greater"}
    
    test[i,"test_ind"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    
    test[i,"test_lot"] <- wilcox.test(par_bac[,paste0(t,"_lot")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    
    
    
    alt <- "less"
    if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,paste0(t,"_lot")])){alt <- "greater"}
    
    test[i,"test_ind_lot"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,paste0(t,"_lot")] , paired = T , alternative = alt)$p.value
    
    
    test[i,"trait"] <- t
    test[i,"p"] <- p
    test[i,"n_ind"] <- nrow(ind)
    test[i,"n_lot"] <- nrow(lot)
    test[i,"nb_geno_ind"] <- length(unique(ind$geno))
    test[i,"nb_geno_lot"] <- length(unique(lot$geno))
    test[i,"nb_geno_alea"] <- length(unique(alea$geno))
    test[i,"i"] <- i_p(p)
    
    i <- i+1
  }
  
  rm(par_bac , ind , lot , alea , tmp)

}



# nb de geno diff?rent entre selection et alea ?
graph <- nb_geno %>% gather(variable , value , c("ind","lot","alea"))
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line()

wilcox.test(nb_geno$lot , nb_geno$alea , alternative = "less")
wilcox.test(nb_geno$ind , nb_geno$alea , alternative = "less")
wilcox.test(nb_geno$lot , nb_geno$ind , alternative = "less")




# hauteur
t <- "hauteur"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)



# N_flag
t <- "N_flag"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# preco
t <- "preco"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# nb_epi
t <- "nb_epi"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# poids_epis
t <- "poids_epis"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)




# en sélectionnant un nombre de grain -------------------------------------

# tri  dans l'ordre décroissant pour que la sélection soit plus facile
pop <- pop %>% arrange(desc(Surface))
moy_geno <- moy_geno %>% arrange(desc(Surface))


test <- data.frame()
don <- data.frame()
nb_geno <- data.frame()

for (n in c(500,400,300,200,100)){
  
  # selection sur lot
  lot <- pop %>% filter(geno == moy_geno$geno[1])
  i <- 1
  while (nrow(lot) < n ){
    i <- i+1
    lot <- rbind(lot , pop %>% filter(geno == moy_geno$geno[i]))
  }
  
  p_lot <- i/nrow(moy_geno)
  i_lot <- i_p(p_lot)
  
  
  # selection sur grain
  ind <- pop[1:n,]
  p_ind <- n/nrow(pop)
  i_ind <- i_p(p_ind)
  
  
  # selection aleatoire
  alea <- pop[sample(row.names(pop) , n),]
  
  
  # selection
  par_bac <- pop %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  
  tmp <- ind %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_ind")] <- tmp[,traits]
  
  tmp <- lot %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_lot")] <- tmp[,traits]
  
  par_bac <- as.data.frame(par_bac)
  
  par_bac$p_lot <- p_lot
  par_bac$i_lot <- i_lot
  par_bac$p_ind <- p_ind
  par_bac$i_ind <- i_ind
  par_bac$n <- n
  
  
  don <- rbind(don,par_bac)
  
  tmp <- data.frame(ind = length(unique(ind$geno)) , lot = length(unique(lot$geno)) , alea = length(unique(alea$geno)) , p_lot = p_lot , i_lot = i_lot , p_ind = p_ind , i_ind = i_ind)
  
  nb_geno <- rbind(nb_geno , tmp)
  
  for (t in c(traits)){
    
    alt <- "less"
    if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,t])){alt <- "greater"}
    
    test[i,"test_ind"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    
    test[i,"test_lot"] <- wilcox.test(par_bac[,paste0(t,"_lot")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    
    
    
    alt <- "less"
    if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,paste0(t,"_lot")])){alt <- "greater"}
    
    test[i,"test_ind_lot"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,paste0(t,"_lot")] , paired = T , alternative = alt)$p.value
    
    
    test[i,"trait"] <- t
    test[i,"p"] <- p
    test[i,"n_ind"] <- nrow(ind)
    test[i,"n_lot"] <- nrow(lot)
    test[i,"nb_geno_ind"] <- length(unique(ind$geno))
    test[i,"nb_geno_lot"] <- length(unique(lot$geno))
    test[i,"nb_geno_alea"] <- length(unique(alea$geno))
    test[i,"i"] <- i_p(p)
    
    i <- i+1
  }
  
  rm(par_bac , ind , lot , alea , tmp)
  
  
}





# nb de geno diff?rent entre selection et alea ?
graph <- nb_geno %>% gather(variable , value , c("ind","lot","alea"))
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line()

wilcox.test(nb_geno$lot , nb_geno$alea , alternative = "less")
wilcox.test(nb_geno$ind , nb_geno$alea , alternative = "less")
wilcox.test(nb_geno$lot , nb_geno$ind , alternative = "less")




# hauteur
t <- "hauteur"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)



# N_flag
t <- "N_flag"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# preco
t <- "preco"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# nb_epi
t <- "nb_epi"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)


# poids_epis
t <- "poids_epis"

graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))

ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations s?lectionn?es" , x = "Intensit? de s?lection") + facet_wrap(~BAC)

test %>% filter(trait == t)









# sélection sur prediction phenomique -------------------------------------

load("resultats_phenomique")

ggplot(resultats_phenomique %>% filter(trait == "hauteur" | trait == "BLUP_hauteur") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats hauteur" , y = "R2" , x = "pré-traitement")

ggplot(resultats_phenomique %>% filter(trait == "preco" | trait == "BLUP_preco") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats preco" , y = "R2" , x = "pré-traitement")
