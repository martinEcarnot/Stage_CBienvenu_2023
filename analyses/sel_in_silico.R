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
traits <- c("preco","N_flag","hauteur","nb_epi","poids_epis","surface_recolte_moy","surface_recolte_min","surface_recolte_max","surface_recolte_moy2","surface_recolte_min2","surface_recolte_max2","PMG","PMG2","poids_moy","poids_min","poids_max","poids_moy2","poids_min2","poids_max2","GSV","GSV2","prot_recolte","nb_grain")

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
        test[i,"NEO"] <- neo
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




# load("../donnees/sel_in_silico")
# 
# sel_in_silico <- rbind(sel_in_silico , test)
# 
# 
# sel_in_silico$NEO <- sel_in_silico$NEO * 177
# 
# 
# save(sel_in_silico , file = "../donnees/sel_in_silico")




load("../donnees/sel_in_silico")

# nb genotypes selectionnes
a <- sel_in_silico %>% filter(trait == "nb_geno") %>% gather(variable , value , c("nb_geno_lot","nb_geno_ind","nb_geno"))

a$NEO <- as.factor(a$NEO)

ggplot(a , aes(x = NEO , y = value , col = variable , fill = variable)) + geom_boxplot() + labs(y = "nb_geno") + facet_wrap(~nsel)


# intensites de selection
a <- sel_in_silico %>% filter(trait == "nb_geno")  %>% gather(variable , value , c("i_ind","i_lot")) %>% mutate(disc = paste(variable , nsel))
ggplot(a , aes(x = NEO , y = value , col = variable)) + geom_point() + geom_line() + labs(y = "intensite de selection" , title = "intensite de selection en fonction du nombre d'epi observe") + facet_wrap(~nsel)











pval <- sel_in_silico %>% group_by(NEO,nsel,trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))


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
 


load("../donnees/sel_in_silico")

neo <- 177
NSEL <- 400


lettres <- sel_in_silico %>% filter(nsel == NSEL & NEO == neo & !trait %in% c("nb_geno",'rdt')) %>% group_by(trait) %>% summarise(pval_ind = mean(t.test_ind, na.rm = T) , pval_lot = mean(t.test_lot , na.rm = T)) %>% gather(variable,value,c("pval_ind","pval_lot"))

lettres$l <- ifelse(lettres$value < 0.001 , "***" , 
                    ifelse(lettres$value < 0.01 , "**",
                           ifelse(lettres$value < 0.05 , "*",NA)))



garde <- c("PMG","GSV2","hauteur","N_flag","nb_epi","poids_epis","preco","prot_recolte","surface_recolte_moy2","surface_recolte_min2","surface_recolte_max","nb_grain")



don <- sel_in_silico %>% 
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
                                                                                                                               ifelse(don$trait == "nb_grain" , "NGE" , ifelse(don$trait == "GSV2" , "GSV" , don$trait)))))))))))

don$trait3 <- factor(don$trait2 , levels = c("PMG","TMG","Taille du plus petit grain","Taille du plus gros grain","PTE","NbEP","NGE","TPG","TPF","PRE","H","GSV"))

for (i in 1:nrow(don)){
  t <- don[i,"trait"]
  v <- don[i,"variable"]
  
  don[i,"htxt"] <- don[i,"conf_haut"] + 0.02 * don[i,"conf_haut"]
}

graph <- don %>% filter(trait3 %in% c("PMG","TMG","PTE","NbEP","NGE","TPG","TPF","PRE","PRE","H","GSV"))
graph$htxt <- graph$conf_haut + 0.02 * graph$conf_haut
graph$fakepoint <- graph$conf_haut + 0.1 * graph$conf_haut
ggplot(graph , aes(x = variable2 , y = value , fill = variable2)) + geom_col() + facet_wrap(~trait3 , scales = "free") + geom_text(aes(label = lettres , y = htxt ) , size = 6) + labs(y = "Progr�s estim�s" , x="" , caption = paste("NEO = ",neo,"et nsel =",NSEL)) + geom_hline(yintercept = 0) + theme(legend.position = "none") + geom_errorbar(aes(ymin = conf_bas , ymax = conf_haut) , width = 0.3) + theme(panel.background = element_blank()) + geom_point(aes(x = variable2 , y = fakepoint) , col = "white" , alpha = 0)





# comparaison avec la theorie ---------------------------------------------


rm(list=setdiff(ls(), c("test2","pop","bac","moy_geno","ngl","n_sel")))



load("../donnees/sel_in_silico")

don <- sel_in_silico %>% filter(trait == "surface_recolte_moy") %>% rename(R_rela_lot = "Rp100_lot" , R_rela_ind = "Rp100_ind") %>% mutate(NEO = as.factor(NEO*177)) 

hline <- don %>% group_by(nsel) %>% summarise(R_ind = mean(R_rela_ind))

ggplot(don , aes(x = NEO , y = R_rela_lot)) + geom_boxplot() + geom_hline(data = hline , aes(yintercept = R_ind , col = "R% grain")) + facet_wrap(~nsel) + labs(x = "Nombre d'épi observé" , y = "progrès génétique en % de l'écart-type" , title = "Comparaison de la sélection sur grain et de la sélection sur lot pour la taille du grain") + theme(axis.text.x = element_text(angle = 90 , size = 7) , legend.title = element_blank())







don <- sel_in_silico %>% filter(trait == "surface_recolte_moy") %>% mutate(RR_exp = R_lot/R_ind) %>% group_by(nsel,NEO) %>% summarise(RR_exp = mean(RR_exp) , Rlot_exp = mean(R_lot) , Rind_exp = mean(R_ind)) %>% as.data.frame()

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


n_sel <- unique(sel_in_silico$nsel)
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


