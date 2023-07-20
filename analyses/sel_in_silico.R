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

moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

rm(opto)

pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F & appel == "present" & semis == "06/01")

pop$pop <- "non_sel"

# traits dont on regarde l'Ã©volution

ggplot(pop , aes(sample = preco)) + geom_qq() + geom_qq_line(col = "red")
ggplot(popsel , aes(sample = preco)) + geom_qq() + geom_qq_line(col = "red")

ggplot(pop , aes(sample = hauteur)) + geom_qq() + geom_qq_line(col = "red")
ggplot(popsel , aes(sample = hauteur)) + geom_qq() + geom_qq_line(col = "red")

ggplot(pop , aes(sample = N_flag)) + geom_qq() + geom_qq_line(col = "red")
ggplot(popsel , aes(sample = N_flag)) + geom_qq() + geom_qq_line(col = "red")

ggplot(pop , aes(sample = nb_epi)) + geom_qq() + geom_qq_line(col = "red")
ggplot(popsel , aes(sample = nb_epi)) + geom_qq() + geom_qq_line(col = "red")

ggplot(pop , aes(sample = poids_epis)) + geom_qq() + geom_qq_line(col = "red")
ggplot(popsel , aes(sample = poids_epis)) + geom_qq() + geom_qq_line(col = "red")




traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis")


# SÃ©lection light pour bonne estimation -----------------------------------

# sur grain individuel

res <- data.frame()
i <- 1

for (p in c(0.5,0.45,0.4,0.35,0.3,0.25,0.2)){
  
  a <- quantile(pop$Surface , probs = 1-p)
  ind <- pop %>% filter(Surface > a)
  
  gen <- moy_geno[which(moy_geno$Surface > quantile(moy_geno$Surface , probs = 1-p)) , "geno"]
  lot <- pop %>% filter(geno %in% gen$geno)
  rm(gen)
  
  alea <- pop[sample(nrow(pop) , round(p*nrow(pop))),]
  
  par_bac <- pop %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  
  tmp <- ind %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_ind")] <- tmp[,traits]
  
  tmp <- lot %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
  par_bac[,paste0(traits,"_lot")] <- tmp[,traits]
  
  par_bac <- as.data.frame(par_bac)
  
  
  
  
  for (t in traits){
    res[i,"sans selection"] <- mean(pop[,t] , na.rm = T)
    res[i,"selection sur grain"] <- mean(ind[,t] , na.rm = T)
    res[i,"selection sur lot"] <- mean(lot[,t] , na.rm = T)
    res[i,"selection aleatoire"] <- mean(alea[,t] , na.rm = T)
    
    res[i,"sd_pop"] <- sd(pop[,t] , na.rm = T)
    res[i,"sd_ind"] <- sd(ind[,t] , na.rm = T)
    res[i,"sd_lot"] <- sd(lot[,t] , na.rm = T)
    res[i,"sd_alea"] <- sd(alea[,t] , na.rm = T)
    
    R <- res[i,"selection sur grain"] - res[i,"sans selection"]
    
    alt <- "less"
    if (R > 0){alt <- "greater"}
    
    # entre sÃ©lections et population de base
    res[i,"t.test_ind"] <- t.test(ind[,t] , pop[,t] , alternative = alt)$p.value
    res[i,"t.test_lot"] <- t.test(lot[,t] , pop[,t] , alternative = alt)$p.value
    
    res[i,"w.test_ind"] <- wilcox.test(ind[,t] , pop[,t] , alternative = alt)$p.value
    res[i,"w.test_lot"] <- wilcox.test(lot[,t] , pop[,t] , alternative = alt)$p.value
    
    res[i,"paired.w.test_ind"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    res[i,"paired.w.test_lot"] <- wilcox.test(par_bac[,paste0(t,"_lot")] , par_bac[,t] , paired = T , alternative = alt)$p.value
    
    
    # comparaison des methodes de selection (sert à rien en vrai)
    # alt <- "less"
    # if (res[i,"R_ind"] > res[i,"R_lot"]){alt <- "greater"}
    # 
    # res[i,"t.test_comp"] <- t.test(ind[,t] , lot[,t] , alternative = alt)$p.value
    # 
    # res[i,"paired.t.test_comp"] <- t.test(par_bac[,paste0(t,"_ind")] , par_bac[,paste0(t,"_lot")] , paired = T , alternative = alt)$p.value
    # 
    # res[i,"w.test_comp"] <- wilcox.test(ind[,t] , lot[,t] , alternative = alt)$p.value
    # 
    # res[i,"paired.w.test_comp"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,paste0(t,"_lot")] , paired = T , alternative = alt)$p.value
    

    res[i,"trait"] <- t
    res[i,"p"] <- p
    res[i,"n_ind"] <- nrow(ind)
    res[i,"n_lot"] <- nrow(lot)
    res[i,"nb_geno_ind"] <- length(unique(ind$geno))
    res[i,"nb_geno_lot"] <- length(unique(lot$geno))
    res[i,"nb_geno_alea"] <- length(unique(alea$geno))
    res[i,"i"] <- i_p(p)
    i <- i+1
  }

}



# nb de geno différent entre selection et alea ?
wilcox.test(unique(res$nb_geno_lot) , unique(res$nb_geno_alea) , alternative = "less")
wilcox.test(unique(res$nb_geno_ind) , unique(res$nb_geno_alea) , alternative = "less")

#wilcox.test(unique(res$nb_geno_lot) , unique(res$nb_geno_alea) , alternative = "less" , paired = T)
#wilcox.test(unique(res$nb_geno_ind), unique(res$nb_geno_alea) , alternative = "less" , paired = T)




# Nombre de genotypes selectionnes
graph <- res %>% filter(trait == "hauteur") %>% gather(variable , value , 19:21) # on prend un trait psq sinon y'a plein de fois la même info
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Nombre de génotypes dans les populations" , x = "Intensité de sélection")


# hauteur
graph <- res %>% filter(trait == "hauteur") %>% gather(variable , value , 1:4)
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations sélectionnées" , x = "Intensité de sélection")


# N_flag
graph <- res %>% filter(trait == "N_flag") %>% gather(variable , value , 1:4)
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations sélectionnées" , x = "Intensité de sélection")


# preco
graph <- res %>% filter(trait == "preco") %>% gather(variable , value , 1:4)
ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations sélectionnées" , x = "Intensité de sélection")











p <- 0.35

a <- quantile(pop$Surface , probs = 1-p)
ind <- pop %>% filter(Surface > a) %>% mutate(pop = "sel_grain")

gen <- moy_geno[which(moy_geno$Surface > quantile(moy_geno$Surface , probs = 1-p)) , "geno"] 

lot <- pop %>% filter(geno %in% gen$geno) %>% mutate(pop = "sel_lot")

rm(gen)

don <- rbind(ind,lot,pop)
don$pop <- as.factor(don$pop)

mod <- lm(hauteur ~ pop , data = don)

lettres <- as.data.frame(cld(glht(model = mod , mcp(pop = "Tukey")) , level = 0.05)$mcletters$Letters) %>% rownames_to_column(var = "pop")

names(lettres)[2] <- "groups"


ggplot(don, aes(x=pop, y=hauteur , col=pop , fill=pop))+
  geom_boxplot(alpha = 0.5)+
  theme(legend.position="none") +
  geom_text(data = lettres, aes(label = groups, y = 150 ), size=5) +
  labs(title = "Hauteur de la pop" , x = "Population" , y = "Hauteur")


t.test(lot$hauteur,pop$hauteur , alternative = "greater")
t.test(ind$hauteur,pop$hauteur , alternative = "greater")

# REGARDER AVEC SCRIPT TEST DE TUCKEY


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















# en sÃ©lectionnant un nombre de grain -------------------------------------

# tri  dans l'ordre dÃ©croissant pour que la sÃ©lection soit plus facile
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











# sÃ©lection sur prediction phenomique -------------------------------------

load("resultats_phenomique")

ggplot(resultats_phenomique %>% filter(trait == "hauteur" | trait == "BLUP_hauteur") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats hauteur" , y = "R2" , x = "prÃ©-traitement")

ggplot(resultats_phenomique %>% filter(trait == "preco" | trait == "BLUP_preco") , aes(x=pretraitement , y = accuracy^2)) + geom_boxplot() + facet_grid(model~donnees) + labs(title = "Resultats preco" , y = "R2" , x = "prÃ©-traitement")
