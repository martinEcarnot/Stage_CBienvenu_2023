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

pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F & appel == "present" & semis == "06/01") %>% mutate(selection = "NON") %>% arrange(desc(Surface))

moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

moy_geno <- moy_geno[which(moy_geno$geno %in% unique(pop$geno)),]

rm(opto)




# traits dont on regarde l'évolution --------------------------------------


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




traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis","surface_recolte_moy","surface_recolte_min","surface_recolte_max","PMG","PMG2","poids_moy","poids_min","poids_max")


# Sélection light pour bonne estimation -----------------------------------

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


resultats <- function(t){
  graph <- don %>% select(t, paste0(t,"_ind") , paste0(t,"_lot") , BAC , i) %>% gather(variable , value , c(t, paste0(t,"_ind") , paste0(t,"_lot"))) %>% mutate(group = paste0(BAC,variable))
  
  g <- ggplot(graph , aes(x = i , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "Moyenne des populations selectionnees" , x = "Intensite de selection" , title = t) + facet_wrap(~BAC)
  
  print(g)
  
  
  print(test %>% filter(trait == t))
  
  f <- as.formula(paste(t,"~ Surface + BAC + geno"))
  
  mod <- lm(f , data = pop)
  drop1(mod , .~. , test = "F")
  
}






# hauteur
resultats("hauteur")



# N_flag
resultats("N_flag")


# preco
resultats("preco")


# nb_epi
resultats("nb_epi")


# poids_epis
resultats("poids_epis")



# surface_recolte_moy
resultats("surface_recolte_moy")


# surface_recolte_max
resultats("surface_recolte_max")


# surface_recolte_min
resultats("surface_recolte_min")



# PMG
resultats("PMG")


# PMG2
resultats("PMG2")


# poids_min
resultats("poids_min")


# poids_max
resultats("poids_max")




# En selectionnant comme dans eq theorique ---------------------------------


mean(table(pop$geno))
hist(table(pop$geno))

# en moyenne, dans la pop non selectionnee un lot a ete mis en terre 5 fois
# on prend ça comme taux de multi du coup

i <- 1
test <- data.frame()


# on fait par bac a chaque fois

ngl <- 5 # nombre de grains par lot, equivalent a nb de grain par epi dans eq theorique
traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis","surface_recolte_moy","surface_recolte_min","surface_recolte_max","surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","PMG","PMG2","poids_moy","poids_min","poids_max","poids_moy2","poids_min2","poids_max2","GSV","GSV2","prot_recolte")

n_sel <- seq(100,600,50)




for (r in 1:100) {

  for (nsel in n_sel){
    
    # caracteristiques de la population selectionnee sur grains
    sel_ind <- pop[1:nsel,]
    
    sel_ind$selection <- "IND"
    
    # intensite de selection ind
    i_ind <- i_p(nsel/nrow(pop))
    
    
    par_bac <- pop %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
    
    tmp <- sel_ind %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
    par_bac[,paste0(traits,"_ind")] <- tmp[,traits]
    
    tmp <- pop %>% group_by(BAC) %>% summarise(rdt = sum(poids_epis , na.rm = T))
    par_bac[,"rdt"] <- tmp[,"rdt"]
    
    tmp <- sel_ind %>% group_by(BAC) %>% summarise(rdt = sum(poids_epis , na.rm = T))
    par_bac[,"rdt_ind"] <- tmp[,"rdt"]
    
    a <- pop %>% group_by(BAC,geno) %>% summarise(count = 1) %>% group_by(BAC) %>% summarise(nb_geno = sum(count)) %>% merge(par_bac , by = "BAC")
    
    b <- sel_ind %>% group_by(BAC,geno) %>% summarise(count = 1) %>% group_by(BAC) %>% summarise(nb_geno_ind = sum(count)) %>% merge(a , by = "BAC")
    
    par_bac <- b
    rm(a,b)
    
    
    for (neo in c(seq(round(nsel/ngl) , nrow(moy_geno) , 10) , nrow(moy_geno))){
      
      # choix au hasard des lots observés
      obs <- sample(row.names(moy_geno) , neo)
      
      # constitution de la population de lots observes
      pop_lot <- moy_geno[obs,] %>% arrange(desc(Surface)) #population d'epi observe sur laquelle on selectionne
      
      # sélection des meilleurs lots dans la population de lots observes
      sel_pop_lot <- pop_lot[1:round(nsel/ngl),]
      
      # caracteristiques de la population selectionnee par lot
      sel_lot <- pop[which(pop$geno %in% sel_pop_lot$geno),]
      
      # intensite de selection lot
      i_lot <- i_p(nrow(sel_pop_lot)/neo)
      
      
      # donnees pour test de wilcoxon apparié
      ## ajout des donnes par bac
      tmp <- sel_lot %>% group_by(BAC) %>% summarise_at(.vars = traits , .funs = mean , na.rm = T)
      par_bac[,paste0(traits,"_lot")] <- tmp[,traits]
      a <- as.data.frame(par_bac)
      
      a$nb_geno_lot <- NULL
      par_bac <- sel_lot %>% group_by(BAC,geno) %>% summarise(count = 1) %>% group_by(BAC) %>% summarise(nb_geno_lot = sum(count)) %>% merge(a , by = "BAC")
      a <- as.data.frame(par_bac)
      
      a$rdt_lot <- NULL
      par_bac <- sel_lot %>% group_by(BAC) %>% summarise(rdt_lot = sum(poids_epis , na.rm = T)) %>% merge(a , by = "BAC")
      
      
      rm(a)
      
      
      
      # donnees pour test student
      sel_lot$selection <- "LOT"
      don_t <- rbind(pop,sel_ind,sel_lot)
      don_t$selection <- relevel(as.factor(don_t$selection) , ref = "NON")
      
      
      # tests et remplissage du tableau de resultats
      
      for (t in c(traits,"nb_geno","rdt")){
        
        test[i,"trait"] <- t
        test[i,"i_ind"] <- i_ind
        test[i,"i_lot"] <- i_lot
        test[i,"NEO"] <- neo/nrow(moy_geno)
        test[i,"nsel"] <- nsel
        test[i,"rep"] <- r
        
        # w.test entre selection et pas de selection
        alt <- "less"
        if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,t])){alt <- "greater"}
        
        test[i,"w.test_ind"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,t] , paired = T , alternative = alt)$p.value
        
        test[i,"w.test_lot"] <- wilcox.test(par_bac[,paste0(t,"_lot")] , par_bac[,t] , paired = T , alternative = alt)$p.value
        
        
        # w.test entre selection sur lot et selection sur grain
        alt <- "less"
        if (mean(par_bac[,paste0(t,"_ind")]) > mean(par_bac[,paste0(t,"_lot")])){alt <- "greater"}
        
        test[i,"w.test_ind_lot"] <- wilcox.test(par_bac[,paste0(t,"_ind")] , par_bac[,paste0(t,"_lot")] , paired = T , alternative = alt)$p.value
      
      
      
      # t.test
      if (t != "nb_geno" & t != "rdt"){
        f <- as.formula(paste(t,"~ BAC + selection"))
        model <- lm(f , data = don_t)
        mod <- summary(model)$coefficients
        conf <- confint(model)
        
        test[i,"t.test_ind"] <- mod["selectionIND","Pr(>|t|)"]
        test[i,"t.test_lot"] <- mod["selectionLOT","Pr(>|t|)"]
        test[i,"R_ind"] <- mod["selectionIND","Estimate"]
        test[i,"R_lot"] <- mod["selectionLOT","Estimate"]
        test[i,"conf_ind_bas"] <- conf["selectionIND",1]
        test[i,"conf_ind_haut"] <- conf["selectionIND",2]
        test[i,"conf_lot_bas"] <- conf["selectionLOT",1]
        test[i,"conf_lot_haut"] <- conf["selectionLOT",2]
        test[i,"Rp100_ind"] <- mod["selectionIND","Estimate"] * 100 / sd(pop[,t] , na.rm = T)
        test[i,"Rp100_lot"] <- mod["selectionLOT","Estimate"] * 100 / sd(pop[,t] , na.rm = T)
        test[i,"confp100_ind_bas"] <- conf["selectionIND",1] * 100 / sd(pop[,t] , na.rm = T)
        test[i,"confp100_ind_haut"] <- conf["selectionIND",2] * 100 / sd(pop[,t] , na.rm = T)
        test[i,"confp100_lot_bas"] <- conf["selectionLOT",1] * 100 / sd(pop[,t] , na.rm = T)
        test[i,"confp100_lot_haut"] <- conf["selectionLOT",2] * 100 / sd(pop[,t] , na.rm = T)
      
        }
        
        
        if (t == "nb_geno"){
          test[i,"nb_geno"] <- mean(par_bac$nb_geno)
          test[i,"nb_geno_ind"] <- mean(par_bac$nb_geno_ind)
          test[i,"nb_geno_lot"] <- mean(par_bac$nb_geno_lot)
        }
        
        if (t == "rdt"){
          test[i,"rdt"] <- mean(par_bac$rdt)
          test[i,"rdt_ind"] <- mean(par_bac$rdt_ind)
          test[i,"rdt_lot"] <- mean(par_bac$rdt_lot)
        }
        
        
        
        
        
        
        
        i <- i+1
      
        
      
      }
      
    }
  
  }
  print(r)
}



sel_in_silico <- test


#save(sel_in_silico , file = "../donnees/sel_in_silico")

test2 <- test
test <- test2 %>% filter(nsel == 400)

# nb genotypes selectionnes
a <- sel_in_silico %>% filter(trait == "nb_geno") %>% gather(variable , value , c("nb_geno_lot","nb_geno_ind","nb_geno"))

a$NEO <- as.factor(a$NEO)

ggplot(a , aes(x = NEO , y = value , col = variable , fill = variable)) + geom_boxplot() + labs(y = "nb_geno") + facet_wrap(~nsel)


# intensites de selection
a <- test %>% filter(trait == "nb_geno")  %>% gather(variable , value , c("i_ind","i_lot"))
ggplot(a , aes(x = NEO , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "intensite de selection" , title = "intensite de selection en fonction du nombre d'epi observe")








resultat <- function(t , NSEL , neo){
  a <- sel_in_silico %>% filter(trait == t & nsel == NSEL & NEO == neo) #%>% gather(variable,value,c("R%_ind","R%_lot"))
  
  b <- ggplot(a , aes(x = i_lot , y = value , col = variable))  + geom_boxplot() + labs(x = "intensite de selection sur epi" , y = t , title = paste("progrès effectué pour",t,"en % de l'ecart type")) + facet_wrap(~nsel) 
  
  print(b)
  
  print(sel_in_silico %>% filter(trait == t & nsel == nsel) %>% select(!c("i_ind","i_lot","nb_geno","nb_geno_ind","nb_geno_lot")) %>% filter(w.test_ind < 0.05))
}


sel_in_silico$NEO <- sel_in_silico$NEO * 177


resultat(t="surface_recolte_moy" , sel_in_silico) #

resultat(t="PMG") #

resultat(t="PMG2") #

resultat(t="hauteur") #?

resultat(t="preco")

resultat(t="N_flag")

resultat(t="nb_epi")

resultat(t="poids_epis") #?

resultat(t="surface_recolte_min") #

resultat(t="surface_recolte_max") #

resultat(t="poids_moy")

resultat(t="poids_min")

resultat(t="poids_max")






pval <- sel_in_silico %>% group_by(NEO,nsel,trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))


signif <- function(t){
don <- pval %>% filter(trait == t)
ggplot(don , aes(x = NEO , y = value , col = variable)) + geom_point() + geom_line() + facet_wrap(~nsel) + geom_hline(yintercept = 0.05 , col = "red")
}






signif(t="surface_recolte_moy") #

signif(t="PMG") #

signif(t="PMG2") #

signif(t="hauteur") #?

signif(t="preco")

signif(t="N_flag")

signif(t="nb_epi")

signif(t="poids_epis") #?

signif(t="surface_recolte_min") #

signif(t="surface_recolte_max") #

signif(t="poids_moy")

signif(t="poids_min")

signif(t="poids_max")



lettres <- sel_in_silico %>% filter(nsel == 400 & NEO == 177 & trait != "nb_geno") %>% group_by(trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))

lettres$l <- ifelse(lettres$value < 0.05 , "*" , NA)


don <- sel_in_silico %>% filter(nsel == 400 & NEO == 100 & trait != "nb_geno") %>% group_by(trait) %>% summarise(R_ind = mean(`R%_ind`, na.rm = T) , R_lot = mean(`R%_lot` , na.rm = T)) %>% gather(variable,value,c("R_ind","R_lot"))

don$lettres <- lettres$l


ggplot(don , aes(x = variable , y = value , fill = variable)) + geom_col() + facet_wrap(~trait) + geom_text(aes(label = lettres , y = value + 1  ))



# comparaison avec la theorie ---------------------------------------------


rm(list=setdiff(ls(), c("test2","pop","bac","moy_geno","ngl","n_sel")))



load("../donnees/sel_in_silico")

don <- sel_in_silico %>% filter(trait == "surface_recolte_moy") %>% rename(R_rela_lot = "R%_lot" , R_rela_ind = "R%_ind") %>% mutate(NEO = as.factor(NEO*177)) 

hline <- don %>% group_by(nsel) %>% summarise(R_ind = mean(R_rela_ind))

ggplot(don , aes(x = NEO , y = R_rela_lot)) + geom_boxplot() + geom_hline(data = hline , aes(yintercept = R_ind , col = "R% grain")) + facet_wrap(~nsel) + labs(x = "Nombre d'épi observé" , y = "progrès génétique en % de l'écart-type" , title = "Comparaison de la sélection sur grain et de la sélection sur lot pour la taille du grain") + theme(axis.text.x = element_text(angle = 90 , size = 7) , legend.title = element_blank())







don <- sel_in_silico %>% filter(trait == "surface_recolte_moy") %>% mutate(RR_exp = R_lot/R_ind , NEO = NEO*177) %>% group_by(nsel,NEO) %>% summarise(RR_exp = mean(RR_exp) , Rlot_exp = mean(R_lot) , Rind_exp = mean(R_ind)) %>% as.data.frame()

#don <- sel_in_silico %>% filter(trait == "surface_recolte_moy") %>% group_by(nsel,NEO) %>% summarise(R_ind = mean(R_ind) , R_lot = mean(R_lot)) %>% mutate(RR_exp = R_lot/R_ind , NEO = NEO*177) %>% as.data.frame()

row.names(don) <- paste(don$NEO,don$nsel)




# Rgrain/Repi
RR <- function(NEO , NGO , vg , vinter , vintra , vpos , NGE , nsel){
  (NGE * NEO * sqrt(vg+vinter+vpos+vintra) * exp((qnorm(1-nsel/NGO)^2 - qnorm(1-nsel/(NGE*NEO))^2) / 2)) / (NGO * sqrt(vg+vinter+vpos+vintra/NGE))
}

Rg <- function(NGO , vg , vinter , vintra , vpos , nsel){
  vg * exp(-(qnorm(1-nsel/NGO)^2)/2) * NGO / (nsel * sqrt(2*pi) * sqrt(vg+vinter+vpos+vintra))
}

Re <- function(NEO , vg , vinter , vintra , vpos , NGE , nsel){
  vg * exp(-(qnorm(1-nsel/(NGE*NEO))^2)/2) * NGE * NEO / (nsel* sqrt(2*pi) * sqrt(vg+vinter+vpos+vintra/NGE))
}



NG_O <- nrow(pop)

# variances estimees dans estimation_variance.R pour la taille des grains

vg <- 1.303^2
vinter <- 1.423^2
vpos <- 1.229^2
vintra <- 2.677^2

ngl <- 5


n_sel <- seq(100,600,50)
RR_test <- data.frame()
i <- 1

for (nsel in n_sel){
  for (neo in c(seq(round(nsel/ngl) , nrow(moy_geno) , 10) , nrow(moy_geno))){
    
    RR_test[i,"RR_th"] <- RR(NEO=neo , NGO=NG_O , vg=vg , vinter=vinter , vintra=vintra , vpos=vpos , NGE=ngl , nsel=nsel)
    
    RR_test[i,"Rlot_th"] <- Re(NEO=neo , vg=vg , vinter=vinter , vintra=vintra , vpos=vpos , NGE=ngl , nsel=nsel)
    
    RR_test[i,"Rind_th"] <- Rg(NGO=NG_O , vg=vg , vinter=vinter , vintra=vintra , vpos=vpos , nsel=nsel)
    
    RR_test[i,"NEO"] <- neo
    
    RR_test[i,"nsel"] <- nsel
    
    i <- i+1
  }
}


row.names(RR_test) <- paste(RR_test$NEO,RR_test$nsel)
RR_test$nsel <- RR_test$NEO <- NULL

RR_test <- merge(RR_test, don , by = "row.names") %>% mutate(nsel = as.factor(nsel))

{
RR_test$aby <- RR_test$abx <- seq(min(RR_test$RR_th) , max(RR_test$RR_th) , length = nrow(RR_test))
ggplot(RR_test , aes(x = RR_th , y = RR_exp)) + geom_point() + geom_smooth(method = "lm" , se = F , aes(col = "Regression")) + geom_line(aes(x = abx , y = aby , col = "y=x") , linewidth = 1) + labs(x = "Repi/Rgrain theorique" , y = "Repi/Rgrain empirique" , title = "Comparaison entre l'approche théorique et la sélection in silico") + scale_colour_manual(values = c( "#619CFF","#F8766D"))
}

{
RR_test$aby <- RR_test$abx <- seq(min(RR_test$Rlot_th) , max(RR_test$Rlot_th) , length = nrow(RR_test))
ggplot(RR_test , aes(x = Rlot_th , y = Rlot_exp)) + geom_point() + geom_smooth(method = "lm" , se = F , aes(col = "Regression")) + geom_line(aes(x = abx , y = aby , col = "y=x") , linewidth = 1) + labs(x = "Repi theorique" , y = "Repi empirique" , title = "Comparaison entre l'approche théorique et la sélection in silico") + scale_colour_manual(values = c( "#619CFF","#F8766D"))
}

{
RR_test$aby <- RR_test$abx <- seq(min(RR_test$Rind_th) , max(RR_test$Rind_th) , length = nrow(RR_test))
ggplot(RR_test , aes(x = Rind_th , y = Rind_exp)) + geom_point() + geom_smooth(method = "lm" , se = F , aes(col = "Regression")) + geom_line(aes(x = abx , y = aby , col = "y=x") , linewidth = 1) + labs(x = "Rgrain theorique" , y = "Rgrain empirique" , title = "Comparaison entre l'approche théorique et la sélection in silico") + scale_colour_manual(values = c( "#619CFF","#F8766D"))
}

cor(RR_test$RR_th , RR_test$RR_exp)
cor(RR_test$RR_th , RR_test$RR_exp)^2

mod <- lm(RR_th~RR_exp , data = RR_test)
summary(mod)
