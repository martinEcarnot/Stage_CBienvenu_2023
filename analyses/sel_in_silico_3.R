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

ngl <- 5 
n_sel <- seq(100,600,50)




for (r in 1:100) {
  
  for (nsel in n_sel){
    
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
      don_t <- rbind(pop,sel_ind,sel_lot)
      don_t$selection <- relevel(as.factor(don_t$selection) , ref = "NON")
    
      
      test[i,"i_ind"] <- i_ind
      test[i,"i_lot"] <- i_lot
      test[i,"NEO"] <- neo
      test[i,"nsel"] <- nsel
      test[i,"rep"] <- r
      
      
      
      # t.test
      model <- lm(surface_recolte_moy ~ selection , data = don_t)
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
      test[i,"Rp100_ind"] <- mod["selectionIND","Estimate"] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      test[i,"Rp100_lot"] <- mod["selectionLOT","Estimate"] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      test[i,"confp100_ind_bas"] <- conf["selectionIND",1] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      test[i,"confp100_ind_haut"] <- conf["selectionIND",2] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      test[i,"confp100_lot_bas"] <- conf["selectionLOT",1] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      test[i,"confp100_lot_haut"] <- conf["selectionLOT",2] * 100 / sd(pop[,"surface_recolte_moy"] , na.rm = T)
      
      i <- i+1
      
      
        
      
      
    }
    
  }
  print(r)
}



# load("../donnees/sel_in_silico")
# 
# sel_in_silico <- rbind(sel_in_silico , test)
# 
# 
# sel_in_silico$NEO <- sel_in_silico$NEO * 177
# 
# 
# save(sel_in_silico , file = "../donnees/sel_in_silico")





# comparaison avec la theorie ---------------------------------------------

don <- test %>% mutate(RR_exp = R_lot/R_ind) %>% group_by(nsel,NEO) %>% summarise(RR_exp = mean(RR_exp) , Rlot_exp = mean(R_lot) , Rind_exp = mean(R_ind)) %>% as.data.frame()

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


n_sel <- unique(don$nsel)
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


