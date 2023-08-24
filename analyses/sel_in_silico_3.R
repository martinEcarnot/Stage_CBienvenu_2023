rm(list=ls())

setwd("~/Stage/Analyses")





library(tidyverse)
library(ggplot2)
#library(multcomp)

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


load("../donnees/opto_recolte_bac")
pop_non_sel_ind <- opto_recolte_bac[,c("Surface","BAC","ind")]
pop_non_sel_ind$selection <- "NON"
rm(opto_recolte_bac)

load("../donnees/bac")
pop <- bac %>% filter(geno != "INCONNU" & is.na(Surface)==F & appel == "present" & semis == "06/01") %>% mutate(selection = "NON") %>% arrange(desc(Surface)) 
pop_non_sel_lot <- pop %>% select(BAC , surface_recolte_moy , selection , geno) %>% rename(surface_moy = surface_recolte_moy)

load("../donnees/opto")
moy_geno <- opto %>% group_by(geno) %>% summarise(Surface = mean(Surface , na.rm = T))
moy_geno <- moy_geno[which(moy_geno$geno %in% unique(pop$geno)),] %>% rename(surface_moy = Surface)
rm(opto)





# En selectionnant comme dans eq theorique ---------------------------------


mean(table(pop$geno))

# en moyenne, dans la pop non selectionnee un lot a ete mis en terre 5 fois
# on prend ça comme taux de multi du coup

i <- 1
test <- data.frame()


# on fait par bac a chaque fois

ngl <- 5 # nombre de grains par lot, equivalent a nb de grain par epi dans eq theorique


n_sel <- seq(100,600,50)


ecty_lot <- sd(pop_non_sel_lot$surface_moy , na.rm = T)
ecty_ind <- sd(pop_non_sel_ind$Surface , na.rm = T)



for (r in 1:100) {
  
  for (nsel in n_sel){
    
    # caracteristiques de la population selectionnee sur grains
    sel_ind <- row.names(pop)[1:nsel]
    
    pop_sel_ind <- pop_non_sel_ind[which(pop_non_sel_ind$ind %in% sel_ind),]
    
    pop_sel_ind$selection <- "IND"
    
    # intensite de selection ind
    i_ind <- i_p(nsel/nrow(pop))
    
    pop_sel_ind <- rbind(pop_sel_ind , pop_non_sel_ind)
    pop_sel_ind$selection <- relevel(as.factor(pop_sel_ind$selection) , ref = "NON")
    
    mod <- lm(Surface ~ selection + BAC , data = pop_sel_ind)
    conf <- confint(mod)
    a <- summary(mod)$coefficients
    
    
    for (neo in c(seq(round(nsel/ngl) , nrow(moy_geno) , 10) , nrow(moy_geno))){
      
      # choix au hasard des lots observés
      obs <- sample(row.names(moy_geno) , neo)
      
      # constitution de la population de lots observes
      pop_lot <- moy_geno[obs,] %>% arrange(desc(surface_moy)) #population d'epi observe sur laquelle on selectionne
      
      # sélection des meilleurs lots dans la population de lots observes
      sel_lot <- pop_lot[1:round(nsel/ngl),"geno"]$geno
      
      # caracteristiques de la population selectionnee par lot
      pop_sel_lot <- pop_non_sel_lot[which(pop_non_sel_lot$geno %in% sel_lot),]
      
      pop_sel_lot$selection <- "LOT"
      
      pop_sel_lot <- rbind(pop_sel_lot , pop_non_sel_lot)
      
      pop_sel_lot$selection <- relevel(as.factor(pop_sel_lot$selection) , ref = "NON")
      # intensite de selection lot
      i_lot <- i_p(length(sel_lot)/neo)
      
      
      mod <- lm(surface_moy ~ selection + BAC , data = pop_sel_lot)
      conf2 <- confint(mod)
      a2 <- summary(mod)$coefficients
      
      # tests et remplissage du tableau de resultats
        
      test[i,"i_ind"] <- i_ind
      test[i,"i_lot"] <- i_lot
      test[i,"NEO"] <- neo
      test[i,"nsel"] <- nsel
      test[i,"rep"] <- r
      
      
      test[i,"t.test_ind"] <- a["selectionIND","Pr(>|t|)"]
      test[i,"R_ind"] <- a["selectionIND","Estimate"]
      test[i,"Rp100_ind"] <- a["selectionIND","Estimate"] * 100 / ecty_ind
      test[i,"conf_ind_bas"] <- conf["selectionIND",1]
      test[i,"conf_ind_haut"] <- conf["selectionIND",2]
      test[i,"confp100_ind_bas"] <- conf["selectionIND",1] * 100 / ecty_ind
      test[i,"confp100_ind_haut"] <- conf["selectionIND",2] * 100 / ecty_ind
      
      test[i,"t.test_lot"] <- a2["selectionLOT","Pr(>|t|)"]
      test[i,"R_lot"] <- a2["selectionLOT","Estimate"]
      test[i,"conf_lot_bas"] <- conf2["selectionLOT",1]
      test[i,"conf_lot_haut"] <- conf2["selectionLOT",2]
      test[i,"Rp100_lot"] <- a2["selectionLOT","Estimate"] * 100 / ecty_lot
      test[i,"confp100_lot_bas"] <- conf2["selectionLOT",1] * 100 / ecty_lot
      test[i,"confp100_lot_haut"] <- conf2["selectionLOT",2] * 100 / ecty_lot
      
      i <- i+1
        
        
        
      
      
    }
    
  }
  print(r)
}




sel_in_silico3 <- test

save(sel_in_silico3 , file = "../donnees/sel_in_silico3")


# comparaison avec la theorie ---------------------------------------------


rm(list=setdiff(ls(), c("pop","bac","moy_geno","ngl","n_sel")))


load("../donnees/sel_in_silico3")


don <- sel_in_silico3  %>% mutate(RR_exp = R_lot/R_ind) %>% group_by(nsel,NEO) %>% summarise(RR_exp = mean(RR_exp) , Rlot_exp = mean(R_lot) , Rind_exp = mean(R_ind)) %>% as.data.frame()

#don <- sel_in_silico3 %>% filter(trait == "surface_recolte_moy") %>% group_by(nsel,NEO) %>% summarise(R_ind = mean(R_ind) , R_lot = mean(R_lot)) %>% mutate(RR_exp = R_lot/R_ind , NEO = NEO*177) %>% as.data.frame()

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


n_sel <- unique(sel_in_silico3$nsel)
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





rm(list=ls())

neo <- 177
NSEL <- 400

load("../donnees/sel_in_silico3")
graph <- sel_in_silico3 %>% 
  filter(NEO == neo & nsel == NSEL) %>% 
  summarise_all(.funs = mean) %>% 
  gather(variable,value,c("Rp100_ind","Rp100_lot")) %>%
  mutate(conf_bas = ifelse(variable == "Rp100_ind" , confp100_ind_bas , confp100_lot_bas) , conf_haut = ifelse(value == "Rp100_ind" , confp100_ind_haut , confp100_lot_haut)) %>% mutate(lettres = "***")

graph$variable2 <- ifelse(graph$variable == "Rp100_ind" , "R grain" , "R epi")

ggplot(graph , aes(x = variable2 , y = value , fill = variable2)) + geom_col() + geom_text(aes(label = lettres , y = conf_haut + 1 )) + labs(y = "Progrès estimés \n(en % de l'écart-type du trait)" , x="" , title = "Progrès estimé pour la taille du grain \npar sélection in silico" , caption = paste("NEO = ",neo,"et nsel =",NSEL)) + geom_hline(yintercept = 0) + theme(legend.position = "none") + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3)





load("../donnees/sel_in_silico")

sel_in_silico %>% filter(trait == "surface_recolte_moy" & NEO == 177 & nsel == 400)
