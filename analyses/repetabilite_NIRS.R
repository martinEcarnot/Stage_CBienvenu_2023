rm(list = ls())

library(tidyverse)

setwd("C:/Users/bienvenu/Documents/Stage/Analyses")

prot_tmp <- read.table("../donnees/data_brute/Stage_CBienvenu_Prediction_proteines_NIRS_épis.csv" , header = T , sep = ";" , dec = ".") %>% column_to_rownames(var = "Nom") 


l <- grep(pattern = ")" , x = row.names(prot_tmp))


prot2 <- prot_tmp[l,] 

prot2$ind <- sapply(strsplit(row.names(prot2) , "(" , fixed = T) , "[" , 1)




# première mesure de chaque geno
l2 <- setdiff(1:nrow(prot_tmp) , l)

prot3 <- prot_tmp[l2,]
prot3$ind <- row.names(prot3)

# on garde que les geno repetes
prot3 <- prot3 %>% filter(ind %in% prot2$ind)



prot <- rbind(prot2,prot3)

rm(prot2,prot3,prot_tmp)



load("../donnees/estim_var")

prot2 <- filter(prot , !ind %in% estim_var$ind)






mod <- lm(Proteines ~ ind , data = prot2)

mod

summary(mod)

plot(mod)


#R2 a 0.8

anova(mod)
