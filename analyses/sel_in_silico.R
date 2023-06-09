rm(list=ls())

setwd("~/Stage/Analyses")

load("../donnees/bac")

library(tidyverse)

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}



pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F)

ggplot(pop , aes(x = Surface)) + geom_histogram()

# ggplot(pop , aes(sample = Surface)) + geom_qq() + geom_qq_line(col = "red")

load("../donnees/opto")

moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

# ggplot(moy_geno , aes(x = Surface)) + geom_histogram()
# ggplot(moy_geno , aes(sample = Surface)) + geom_qq() + geom_qq_line(col = "red")

# traits dont on regarde l'évolution
traits <- c("preco","N_flag","hauteur")



# en selectionnant selon un seuil -----------------------------------------


# seuils de troncation
seuils <- seq(19, 22 , 0.2)

# Selection sur grain individuel
selection <- data.frame()
j <- 1

for (a in seuils){
  
  # selection sur lot
  gensel <- moy_geno %>% filter(Surface > a)
  popsel_lot <- pop %>% filter(geno %in% gensel$geno)
  pl <- nrow(popsel_lot)/nrow(pop)
  il <- i_p(pl)
  
  
  # selection sur grain
  popsel_grain <- pop %>% filter(Surface > a)
  pg <- nrow(popsel_grain)/nrow(pop)
  ig <- i_p(pg)
  
  
  # selection aleatoire
  alea <- sample(row.names(pop) , nrow(popsel_grain))
  popsel_alea <- pop[alea,]
  
  
  for (t in traits){
    # selection sur lot
    selection[j,"trait"] <- t
    selection[j,"R"] <- mean(popsel_lot[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    selection[j,"nb_geno"] <- length(unique(popsel_lot$geno))
    selection[j,"i"] <- il
    selection[j,"a"] <- a
    
    if (selection[j,"R"] > 0){alt <- "greater"}
    if (selection[j,"R"] < 0){alt <- "less"}
    selection[j,"test_alea"] <- t.test(popsel_lot[,t] , popsel_alea[,t] , alternative = alt)$p.value
    selection[j,"test_pop"] <- t.test(popsel_lot[,t] , pop[,t] , alternative = alt)$p.value
    selection[j,"type"] <- "lot"
    j <- j+1
    
    # selection sur grain
    selection[j,"trait"] <- t
    selection[j,"R"] <- mean(popsel_grain[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    selection[j,"nb_geno"] <- length(unique(popsel_grain$geno))
    selection[j,"i"] <- ig
    selection[j,"a"] <- a
    
    if (selection[j,"R"] > 0){alt <- "greater"}
    if (selection[j,"R"] < 0){alt <- "less"}
    selection[j,"test_alea"] <- t.test(popsel_grain[,t] , popsel_alea[,t] , alternative = alt)$p.value
    selection[j,"test_pop"] <- t.test(popsel_grain[,t] , pop[,t] , alternative = alt)$p.value
    selection[j,"type"] <- "grain"
    j <- j+1
    
    
    # selection aleatoire
    selection[j,"trait"] <- t
    selection[j,"R"] <- mean(popsel_alea[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    selection[j,"nb_geno"] <- length(unique(popsel_alea$geno))
    selection[j,"i"] <- ig
    selection[j,"a"] <- a
    
    if (selection[j,"R"] > 0){alt <- "greater"}
    if (selection[j,"R"] < 0){alt <- "less"}
    selection[j,"test_alea"] <- t.test(popsel_alea[,t] , popsel_alea[,t] , alternative = alt)$p.value
    selection[j,"test_pop"] <- t.test(popsel_alea[,t] , pop[,t] , alternative = alt)$p.value
    selection[j,"type"] <- "alea"
    j <- j+1

  }
}




pour_test <- selection %>% select(trait,R,i,type,nb_geno) %>% filter(type != "lot") %>% reshape(timevar = "type" , idvar = c("i","trait") , direction = "wide" )




ggplot(selection , aes(x = i , y = a , col = type)) + geom_point() + geom_line()

#nb de genotypeq
don <- pour_test %>% filter(trait == "preco") # on prend un trait au hasard c'est pareil pour tous
don2 <- selection %>% filter(trait == "preco")
ggplot(don2 , aes(x = i , y = nb_geno , col = type)) + geom_point() + geom_line() + labs(title = "preco")
t.test(don$nb_geno.grain , don$nb_geno.alea , paired = T , alternative = "less")

ggplot(don , aes(x = nb_geno.grain - nb_geno.alea)) + geom_histogram()
ggplot(don , aes(sample = nb_geno.grain - nb_geno.alea)) + geom_qq() + geom_qq_line(col = "red")

wilcox.test(don$nb_geno.grain , don$nb_geno.alea , paired = T , alternative = "less")

# preco
don <- pour_test %>% filter(trait == "preco") 
don2 <- selection %>% filter(trait == "preco")
ggplot(don2 , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "preco")
t.test(don$R.grain , don$R.alea , paired = T , alternative = "less")

ggplot(don , aes(x = R.grain - R.alea)) + geom_histogram()
ggplot(don , aes(sample = R.grain - R.alea)) + geom_qq() + geom_qq_line(col = "red")

wilcox.test(don$R.grain , don$R.alea , paired = T , alternative = "less")

# N_flag
don <- pour_test %>% filter(trait == "N_flag")
don2 <- selection %>% filter(trait == "N_flag")
ggplot(don2 , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "N_flag") 
t.test(don$R.grain , don$R.alea , paired = T , alternative = "greater")

ggplot(don , aes(x = R.grain - R.alea)) + geom_histogram()
ggplot(don , aes(sample = R.grain - R.alea)) + geom_qq() + geom_qq_line(col = "red")

wilcox.test(don$R.grain , don$R.alea , paired = T , alternative = "greater")

# hauteur
don <- pour_test %>% filter(trait == "hauteur") 
don2 <- selection %>% filter(trait == "hauteur")
ggplot(don2 , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "hauteur")
t.test(don$R.grain , don$R.alea , paired = T , alternative = "greater")

ggplot(don , aes(x = R.grain - R.alea)) + geom_histogram()
ggplot(don , aes(sample = R.grain - R.alea)) + geom_qq() + geom_qq_line(col = "red")

wilcox.test(don$R.grain , don$R.alea , paired = T , alternative = "greater")















# en sélectionnant un nombre de grain -------------------------------------

# tri  dans l'ordre décroissant pour que la sélection soit plus facile
pop <- pop %>% arrange(desc(Surface))
moy_geno <- moy_geno %>% arrange(desc(Surface))

nsel <- seq(200,500,10)


sel_lot <- data.frame()
j <- 1

for (n in nsel){
  
  # selection sur lot
  gensel <- moy_geno[1:round(n/10),]
  popsel_lot <- pop %>% filter(geno %in% gensel$geno)
  il <- i_p(nrow(gensel)/nrow(moy_geno))
  
  
  # selection sur grain
  popsel_grain <- pop[1:n,]
  ig <- i_p(n/nrow(pop))
  
  
  # selection aleatoire
  alea <- sample(row.names(pop) , n)
  popsel_alea <- pop[alea,]
  
  for (t in traits){
    # selection sur lot
    sel_lot[j,"trait"] <- t
    sel_lot[j,"R"] <- mean(popsel_lot[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    sel_lot[j,"nb_geno"] <- length(unique(popsel_lot$geno))
    sel_lot[j,"i"] <- il
    sel_lot[j,"nsel"] <- n
    
    if (sel_lot[j,"R"] > 0){alt <- "greater"}
    if (sel_lot[j,"R"] < 0){alt <- "less"}
    sel_lot[j,"test_alea"] <- t.test(popsel_lot[,t] , popsel_alea[,t] , alternative = alt)$p.value
    sel_lot[j,"test_pop"] <- t.test(popsel_lot[,t] , pop[,t] , alternative = alt)$p.value
    sel_lot[j,"type"] <- "lot"
    j <- j+1
    
    # selection sur grain
    sel_lot[j,"trait"] <- t
    sel_lot[j,"R"] <- mean(popsel_grain[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    sel_lot[j,"nb_geno"] <- length(unique(popsel_grain$geno))
    sel_lot[j,"i"] <- ig
    sel_lot[j,"nsel"] <- n
    
    if (sel_lot[j,"R"] > 0){alt <- "greater"}
    if (sel_lot[j,"R"] < 0){alt <- "less"}
    sel_lot[j,"test_alea"] <- t.test(popsel_grain[,t] , popsel_alea[,t] , alternative = alt)$p.value
    sel_lot[j,"test_pop"] <- t.test(popsel_grain[,t] , pop[,t] , alternative = alt)$p.value
    sel_lot[j,"type"] <- "grain"
    j <- j+1
    
    
    # selection aleatoire
    sel_lot[j,"trait"] <- t
    sel_lot[j,"R"] <- mean(popsel_alea[,t] , na.rm = T) - mean(pop[,t] , na.rm = T)
    sel_lot[j,"nb_geno"] <- length(unique(popsel_alea$geno))
    sel_lot[j,"i"] <- ig
    sel_lot[j,"nsel"] <- n
    
    if (sel_lot[j,"R"] > 0){alt <- "greater"}
    if (sel_lot[j,"R"] < 0){alt <- "less"}
    sel_lot[j,"test_alea"] <- t.test(popsel_alea[,t] , popsel_alea[,t] , alternative = alt)$p.value
    sel_lot[j,"test_pop"] <- t.test(popsel_alea[,t] , pop[,t] , alternative = alt)$p.value
    sel_lot[j,"type"] <- "alea"
    j <- j+1
    
  }
  
}




ggplot(sel_lot , aes(x = i , y = nsel , col = type)) + geom_point() + geom_line()

ggplot(sel_lot , aes(x = i , y = nb_geno , col = type)) + geom_point() + geom_line() + labs(title = "preco")

ggplot(sel_lot %>% filter(trait == "preco") , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "preco")


ggplot(sel_lot %>% filter(trait == "N_flag") , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "N_flag")

ggplot(sel_lot %>% filter(trait == "hauteur") , aes(x = i , y = R , col = type)) + geom_point() + geom_line() + labs(title = "hauteur")






# sélection sur prediction phenomique -------------------------------------

load("resultats_phenomique")

ggplot(resultats_phenomique %>% filter(trait == "hauteur" | trait == "BLUP_hauteur") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats hauteur" , y = "R2" , x = "pré-traitement")

ggplot(resultats_phenomique %>% filter(trait == "preco" | trait == "BLUP_preco") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats preco" , y = "R2" , x = "pré-traitement")
