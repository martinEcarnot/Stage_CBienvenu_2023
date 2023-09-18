rm(list=ls())

setwd("~/Stage/Analyses")


load("../donnees/bac")

library(tidyverse)
library(ggplot2)
#library(multcomp)

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

# en moyenne, dans la pop non selectionnee un lot a ete mis en terre 5 fois
# on prend ça comme taux de multi du coup

i <- 1
test <- data.frame()


# on fait par bac a chaque fois

ngl <- 5 # nombre de grains par lot, equivalent a nb de grain par epi dans eq theorique
traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis","surface_recolte_moy2","PMG2","GSV2","prot_recolte","nb_grain")

n_sel <- seq(100,600,50)




for (r in 1:100) {

  for (nsel in n_sel){
    
    
    #population non sélectionnée
    non_sel <- pop[sample(1:nrow(pop) , nsel) , ]
    
    
    # caracteristiques de la population selectionnee sur grains
    sel_ind <- pop[1:nsel,]
    
    sel_ind$selection <- "IND"
    
    # intensite de selection ind
    i_ind <- i_p(nsel/nrow(pop))
    
    
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
      
      
      
      # donnees pour test student
      sel_lot$selection <- "LOT"
      don_t <- rbind(non_sel,sel_ind,sel_lot)
      don_t$selection <- relevel(as.factor(don_t$selection) , ref = "NON")
      
      
      # tests et remplissage du tableau de resultats
      
      for (t in traits){
        
        test[i,"trait"] <- t
        test[i,"i_ind"] <- i_ind
        test[i,"i_lot"] <- i_lot
        test[i,"NEO"] <- neo
        test[i,"nsel"] <- nsel
        test[i,"rep"] <- r
        
      
      # t.test
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
        
      i <- i+1
      
      }
      
    }
  
  }
  print(r)
}




sel_in_silico_4 <- test

save(sel_in_silico_4 , file = "../donnees/sel_in_silico_4")




load("../donnees/sel_in_silico_4")


pval <- sel_in_silico_4 %>% group_by(NEO,nsel,trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))


signif <- function(t){
don <- pval %>% filter(trait == t)
ggplot(don , aes(x = NEO , y = value , col = variable)) + geom_point() + geom_line() + facet_wrap(~nsel) + geom_hline(yintercept = 0.05 , col = "red")
}


signif(t="surface_recolte_moy") #
signif(t="surface_recolte_moy2") #

signif(t="surface_recolte_min") #
signif(t="surface_recolte_min2")

signif(t="surface_recolte_max") #
signif(t="surface_recolte_max2") #

signif(t="GSV")
signif(t="GSV2") #

signif(t="PMG") #
signif(t="PMG2") #

signif(t="prot_recolte")

signif(t="hauteur") #?

signif(t="preco")

signif(t="N_flag")

signif(t="nb_epi")

signif(t="poids_epis")


signif(t="poids_moy")
signif(t="poids_moy2")

signif(t="poids_min")
signif(t="poids_min2")

signif(t="poids_max")
signif(t="poids_max2")

signif(t="nb_grain")
 


load("../donnees/sel_in_silico_4")

neo <- 177
NSEL <- 250


lettres <- sel_in_silico_4 %>% filter(nsel == NSEL & NEO == neo & !trait %in% c("nb_geno",'rdt')) %>% group_by(trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))

lettres$l <- ifelse(lettres$value < 0.001 , "***" , 
                    ifelse(lettres$value < 0.01 , "**",
                           ifelse(lettres$value < 0.05 , "*",NA)))



garde <- c("PMG2","GSV2","hauteur","N_flag","nb_epi","poids_epis","preco","prot_recolte","surface_recolte_moy2","nb_grain")



don <- sel_in_silico_4 %>% 
  filter(nsel == NSEL & NEO == neo & trait %in% garde) %>% 
  group_by(trait) %>% 
  summarise(R_ind = mean(R_ind, na.rm = T) , R_lot = mean(R_lot , na.rm = T) , conf_ind_bas = mean(conf_ind_bas) , conf_ind_haut = mean(conf_ind_haut) , conf_lot_bas = mean(conf_lot_bas) , conf_lot_haut = mean(conf_lot_haut) , sd_ind = sd(R_ind , na.rm = T) , sd_lot = sd(R_lot , na.rm = T)) %>% 
gather(variable,value,c("R_ind","R_lot")) %>%
  mutate(conf_bas = ifelse(variable == "R_ind" , conf_ind_bas , conf_lot_bas) , conf_haut = ifelse(value == "R_ind" , conf_ind_haut , conf_lot_haut) , conf_sd = ifelse(variable == "R_ind" , sd_ind , sd_lot))

lettres <- lettres %>% filter(trait %in% garde)

don$lettres <- lettres$l

don$variable2 <- ifelse(don$variable == "R_ind" , "R grain" , "R epi")

don$trait2 <- ifelse(don$trait == "surface_recolte_moy2" , "TMG",
                     ifelse(don$trait == "surface_recolte_min2","Taille du plus petit grain",
                            ifelse(don$trait == "surface_recolte_max" , "Taille du plus gros grain",
                                   ifelse(don$trait == "hauteur","H",
                                          ifelse(don$trait == "N_flag","TPF",
                                                 ifelse(don$trait == "prot_recolte","TPG",
                                                        ifelse(don$trait == "nb_epi","NbEP",
                                                               ifelse(don$trait == "poids_epis", "PTE" , ifelse(don$trait == "preco" , "PRE",
                                                                                                                               ifelse(don$trait == "nb_grain" , "NGE" , ifelse(don$trait == "GSV2" , "GSV" , ifelse(don$trait == "PMG2" , "PMG" , don$trait))))))))))))

don$trait3 <- factor(don$trait2 , levels = c("PMG","TMG","Taille du plus petit grain","Taille du plus gros grain","PTE","NbEP","NGE","TPG","TPF","PRE","H","GSV"))

for (i in 1:nrow(don)){
  t <- don[i,"trait"]
  v <- don[i,"variable"]
  
  don[i,"htxt"] <- don[i,"conf_haut"] + 0.02 * don[i,"conf_haut"]
}

graph <- don %>% filter(trait3 %in% c("PMG","TMG","PTE","NbEP","NGE","TPG","TPF","PRE","PRE","H","GSV"))
graph$htxt <- graph$conf_haut + 0.02 * graph$conf_haut
graph$fakepoint <- graph$conf_haut + 0.1 * graph$conf_haut
ggplot(graph , aes(x = variable2 , y = value , fill = variable2)) + geom_col() + facet_wrap(~trait3 , scales = "free") + geom_text(aes(label = lettres , y = htxt ) , size = 6) + labs(y = "Progrès estimés" , x="" , caption = paste("NEO = ",neo,"et nsel =",NSEL)) + geom_hline(yintercept = 0) + theme(legend.position = "none") + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3) + theme(panel.background = element_blank()) + geom_point(aes(x = variable2 , y = fakepoint) , col = "white" , alpha = 0) + scale_fill_manual(values = c("#00BA38" , "#F8766D"))


graph <- don %>% filter(trait3 == "TMG")
graph$fakepoint <- graph$conf_haut + 0.1 * graph$conf_haut
ggplot(graph , aes(x = variable2 , y = value , fill = variable2)) + geom_col() + geom_text(aes(label = lettres , y = htxt ) , size = 6) + labs(y = "Progrès estimés (mm²)" , x="" , caption = paste("NEO = ",neo,"et nsel =",NSEL) , title = "Progrès estimé pour la taille des grains") + geom_hline(yintercept = 0) + theme(legend.position = "none") + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3) + theme(panel.background = element_blank()) + geom_point(aes(x = variable2 , y = fakepoint) , col = "white" , alpha = 0) + scale_fill_manual(values = c("#00BA38" , "#F8766D"))


# comparaison avec la theorie ---------------------------------------------


rm(list=setdiff(ls(), c("test2","pop","bac","moy_geno","ngl","n_sel")))



load("../donnees/sel_in_silico_4")



don <- sel_in_silico_4 %>% filter(trait == "surface_recolte_moy2") %>% mutate(RR_exp = R_lot/R_ind) %>% group_by(nsel,NEO) %>% summarise(RR_exp = mean(RR_exp) , Rlot_exp = mean(R_lot) , Rind_exp = mean(R_ind)) %>% as.data.frame()

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


n_sel <- unique(sel_in_silico_4$nsel)
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

RR_test$RgRe_th <- RR_test$Rind_th / RR_test$Rlot_th
RR_test$RgRe_exp <- RR_test$Rind_exp / RR_test$Rlot_exp



RR_test$one <- ifelse(RR_test$NEO == 177 & RR_test$nsel == 250 , "OUI" , "NON")


{
  mod <- lm(RR_th~RR_exp , data = RR_test)
  
  b <- round(mod$coefficients["(Intercept)"],2)
  a <- round(mod$coefficients["RR_exp"],2)
  r2 <- round(summary(mod)$r.squared,2)
  
ggplot(RR_test , aes(x = RR_exp , y = RR_th)) + geom_point() + geom_smooth(method = "lm" , se = F , col = "#F8766D") + labs(y = "Repi / Rgrain théorique" , x = "Repi / Rgrain empirique") + annotate(geom = "text" , label = paste0("y = ",b," + ",a,"x \nR² = ",r2) , x = 0.3 , y = 1.3 , col = "#F8766D" , size = 5) + annotate("rect", xmin = -0.02 , xmax = 0.6, ymin = 1.1, ymax = 1.5 , alpha = 0 , col = "black") + theme(panel.background = element_blank())
}




{
  mod <- lm(RR_exp~RR_th , data = RR_test)
  
  b <- round(mod$coefficients["(Intercept)"],2)
  a <- round(mod$coefficients["RR_th"],2)
  r2 <- round(summary(mod)$r.squared,2)
  
  ggplot(RR_test , aes(x = RR_th , y = RR_exp ,  col = one , shape = one , size = one)) + geom_point() + geom_smooth(method = "lm" , se = F , col = "#F8766D" , linewidth = 1.2) + labs(y = "Repi / Rgrain empirique" , x = "Repi / Rgrain théorique") + annotate(geom = "text" , label = paste0("y = ",b," + ",a,"x \nR² = ",r2) , x = 0.3 , y = 1.3 , col = "#F8766D" , size = 5) + annotate("rect", xmin = -0.02 , xmax = 0.6, ymin = 1.1, ymax = 1.5 , alpha = 0 , col = "black") + theme(panel.background = element_blank() , legend.position = "none") + scale_colour_manual(values = c("#000000", "green")) + scale_shape_manual(values = c(16,17)) + scale_size_manual(values = c(1.5,4))
}






{
  mod <- lm(RR_th~RR_exp , data = RR_test)
  
  b <- round(mod$coefficients["(Intercept)"],2)
  a <- round(mod$coefficients["RR_exp"],2)
  r2 <- round(summary(mod)$r.squared,2)
  
  ggplot(RR_test , aes(x = RR_exp , y = RR_th)) + geom_point() + geom_smooth(method = "lm" , se = F , col = "#F8766D" , linewidth = 1) + labs(y = "Repi / Rgrain théorique" , x = "Repi / Rgrain empirique") + annotate(geom = "text" , label = paste0("y = ",b," + ",a,"x \nR² = ",r2) , x = 0.3 , y = 1.3 , col = "#F8766D" , size = 5) + annotate("rect", xmin = -0.02 , xmax = 0.6, ymin = 1.1, ymax = 1.5 , alpha = 0 , col = "black") + theme(panel.background = element_blank())

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


