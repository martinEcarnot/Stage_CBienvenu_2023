rm(list=ls())

setwd("~/Stage/Analyses")


load("../donnees/bac")

load("../donnees/opto_recolte_bac")
opto_recolte_bac <- opto_recolte_bac %>% filter(Ref.Ech == "S01") %>% mutate(selection = "IND")

library(tidyverse)

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


load("../donnees/opto")

pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F & appel == "present" & semis == "06/01") %>% mutate(selection = "NON") %>% arrange(desc(Surface))

moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))

moy_geno <- moy_geno[which(moy_geno$geno %in% unique(pop$geno)),]

rm(opto)







# En selectionnant comme dans eq theorique ---------------------------------


mean(table(pop$geno))
hist(table(pop$geno))

# en moyenne, dans la pop non selectionnee un lot a ete mis en terre 5 fois
# on prend ça comme taux de multi du coup

i <- 1
test <- data.frame()


# on fait par bac a chaque fois

ngl <- 5 # nombre de grains par lot, equivalent a nb de grain par epi dans eq theorique


n_sel <- seq(100,600,50)




for (r in 1:100) {

  for (nsel in n_sel){
    
    # caracteristiques de la population selectionnee sur grains
    ind_sel <- row.names(pop[1:nsel,])
    
    pop_sel_ind <- opto_recolte_bac[which(opto_recolte_bac$ind %in% ind_sel),c("selection","Surface","BAC")]
    
    
    # donnees pour calculer R_grain
    calc_rg <- pop %>% select(selection , Surface) %>% mutate(BAC = 7) %>% rbind(pop_sel_ind)
    calc_rg$selection <- relevel(as.factor(calc_rg$selection) , ref = "NON")
    
    mod <- lm(Surface ~ selection + BAC , data = calc_rg)
    
    Rg <- mod$coefficients["selectionIND"]
    t.rg <- summary(mod)$coefficients["selectionIND","Pr(>|t|)"]
    
    
    for (neo in c(seq(round(nsel/ngl) , nrow(moy_geno) , 10) , nrow(moy_geno))){
      
      # choix au hasard des lots observés
      obs <- sample(row.names(moy_geno) , neo)
      
      # constitution de la population de lots observes
      pop_lot <- moy_geno[obs,] %>% arrange(desc(Surface)) #population d'epi observe sur laquelle on selectionne
      
      # sélection des meilleurs lots dans la population de lots observes
      sel_pop_lot <- pop_lot[1:round(nsel/ngl),]
      
      # caracteristiques de la population selectionnee par lot
      sel_lot <- pop[which(pop$geno %in% sel_pop_lot$geno),] %>% select(surface_recolte_moy , BAC) %>% mutate(selection = "LOT") %>% na.omit()
      
      names(sel_lot)[1] <- "Surface"
      
      
      # donnees pour test student
      calc_re <- moy_geno %>% select(Surface) %>% mutate(selection = "NON" , BAC = 7) %>% rbind(sel_lot)
      calc_re$selection <- relevel(as.factor(calc_re$selection) , ref = "NON")
      
      mod <- lm(Surface ~ selection + BAC , data = calc_re)
      
      Re <- mod$coefficients["selectionLOT"]
      t.re <- summary(mod)$coefficients["selectionLOT","Pr(>|t|)"]
      
      # tests et remplissage du tableau de resultats
        
      test[i,"i_ind"] <- i_p(nsel/nrow(pop))
      test[i,"i_lot"] <- i_p(nrow(sel_pop_lot)/neo)
      test[i,"NEO"] <- neo
      test[i,"nsel"] <- nsel
      test[i,"rep"] <- r
      test[i,"Rg"] <- Rg
      test[i,"Re"] <- Re
      test[i,"Rg_p100"] <- Rg / sd(pop$Surface , na.rm = T)
      test[i,"Re_p100"] <- Re / sd(moy_geno$Surface , na.rm = T)
      test[i,"t.rg"] <- t.rg
      test[i,"t.re"] <- t.re
        
      i <- i+1
      
    
      
    }
  
  }
  print(r)
}



sel_in_silico_2 <- test
save(sel_in_silico_2 , file = "../donnees/sel_in_silico_2")








rm(list=ls())

load("../donnees/sel_in_silico_2")

neo <- 100
NSEL <- 400


lettres <- sel_in_silico_2 %>% filter(nsel == NSEL & NEO == neo) %>% gather(variable,value,c("Rg","Re"))

ggplot(lettres , aes(x = variable , y = value)) + geom_boxplot()


don <- sel_in_silico_2 %>% 
  filter(nsel == NSEL & NEO == neo) %>% 
  group_by(trait) %>% 
  summarise(R_ind = mean(Rp100_ind, na.rm = T) , R_lot = mean(Rp100_lot , na.rm = T) , conf_ind_bas = mean(confp100_ind_bas) , conf_ind_haut = mean(confp100_ind_haut) , conf_lot_bas = mean(confp100_lot_bas) , conf_lot_haut = mean(confp100_lot_haut) , sd_ind = sd(Rp100_ind , na.rm = T) , sd_lot = sd(Rp100_lot , na.rm = T)) %>% 
gather(variable,value,c("R_ind","R_lot")) %>%
  mutate(conf_bas = ifelse(variable == "R_ind" , conf_ind_bas , conf_lot_bas) , conf_haut = ifelse(value == "R_ind" , conf_ind_haut , conf_lot_haut) , conf_sd = ifelse(variable == "R_ind" , sd_ind , sd_lot))

lettres <- lettres %>% filter(trait %in% garde)

don$lettres <- lettres$l

don$variable2 <- ifelse(don$variable == "R_ind" , "R grain" , "R epi")

don$trait2 <- ifelse(don$trait == "surface_recolte_moy" , "Taille moyenne des grains",
                     ifelse(don$trait == "surface_recolte_min","Taille du plus petit grain",
                            ifelse(don$trait == "surface_recolte_max" , "Taille du plus gros grain",
                                   ifelse(don$trait == "hauteur","Hauteur",
                                          ifelse(don$trait == "N_flag","Taux N feuille drapeau",
                                                 ifelse(don$trait == "prot_recolte","Taux N grains",
                                                        ifelse(don$trait == "nb_epi","Nombre d'épis",
                                                               ifelse(don$trait == "poids_epis", "Poids total d'épis" , ifelse(don$trait == "preco" , "Précocité",
                                                                                                                               ifelse(don$trait == "nb_grain" , "Nombre de grains par épi" , don$trait))))))))))

don$trait3 <- factor(don$trait2 , levels = c("PMG","Taille moyenne des grains","Taille du plus petit grain","Taille du plus gros grain","Poids total d'épis","Nombre d'épis","Nombre de grains par épi","Taux N grains","Taux N feuille drapeau","Précocité","Hauteur","GSV"))

ggplot(don , aes(x = variable2 , y = value , fill = variable2)) + geom_col() + facet_wrap(~trait3) + geom_text(aes(label = lettres , y = conf_haut + 1 )) + labs(y = "Progrès estimés (en % de l'écart-type du trait)" , x="" , title = "Progrès estimés par sélection in silico" , caption = paste("NEO = ",neo,"et nsel =",NSEL)) + geom_hline(yintercept = 0) + theme(legend.position = "none") + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3)





# comparaison avec la theorie ---------------------------------------------





load("../donnees/sel_in_silico_2")


# Rgrain/Repi
RR <- function(NEO , NGO , vg , vinter , vintra , vpos , NGE , nsel){
  (NGE * NEO * sqrt(vg+vinter+vpos+vintra) * exp((qnorm(1-nsel/NGO)^2 - qnorm(1-nsel/(NGE*NEO))^2) / 2)) / (NGO * sqrt(vg+vinter+vpos+vintra/NGE))
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
    
    RR_test[i,"NEO"] <- neo
    
    RR_test[i,"nsel"] <- nsel
    
    i <- i+1
  }
}


row.names(RR_test) <- paste(RR_test$NEO,RR_test$nsel)
RR_test$nsel <- RR_test$NEO <- NULL

don <- sel_in_silico_2 %>% mutate(RR_exp = Rg/Re) %>% group_by(NEO , nsel) %>% summarise(RR_exp = mean(RR_exp)) %>% as.data.frame()

row.names(don) <- paste(don$NEO,don$nsel)

RR_test <- merge(RR_test, don , by = "row.names") %>% mutate(nsel = as.factor(nsel))



{
  mod <- lm(RR_th~RR_exp , data = RR_test)
  
  b <- round(mod$coefficients["(Intercept)"],2)
  a <- round(mod$coefficients["RR_exp"],2)
  r2 <- round(summary(mod)$r.squared,2)
  
ggplot(RR_test , aes(x = RR_exp , y = RR_th)) + geom_point() + geom_smooth(method = "lm" , se = F , col = "#F8766D") + labs(y = "Repi / Rgrain théorique" , x = "Repi / Rgrain empirique" , title = "Comparaison entre l'approche théorique et la sélection in silico") + annotate(geom = "text" , label = paste0("y = ",b," + ",a,"x \nR² = ",r2) , x = 0.3 , y = 1 , col = "#F8766D" , size = 5) + annotate("rect", xmin = 0.02 , xmax = 0.55, ymin = 0.8, ymax = 1.2 , alpha = 0 , col = "black") + labs()
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


