library(ggplot2)
library(tidyverse)

vg <- 1
vintra <- 2
ve <- 2

NGE <- 70

sepi <- sqrt(vg+ve+vintra/NGE)

p <- 0.05

# sans trade off ----------------------------------------------------------

# population d'épi
epis <- rnorm(mean = 0 , sd = sepi , n = 10000)
seuil <- qnorm(mean = 0 , sd = sepi , p = 1-p)
epis_sel <- epis[which(epis > seuil)]

# population de grains
grains_sel <- c()
for (e in epis_sel){
  grains_sel <- c(grains_sel , rnorm(mean = e , sd = vintra , n = NGE))
}



# avec trade off ----------------------------------------------------------

# pente du trade off
a <- -10
oo <- 70

# population de grains

grains_sel_to <- c()

for (e in epis_sel){
  grains_sel_to <- c(grains_sel_to , rnorm(mean = e , sd = vintra , n = oo + a*e))
}




tab <- data.frame(ech = c(epis , grains_sel , grains_sel_to) , 
                  dist = c(rep("moyenne épis" , length(epis)) , rep("Grains sans trade off" , length(grains_sel)) , rep("Grains avec trade off" , length(grains_sel_to))))

ggplot(tab , aes(x = ech , col = dist , fill = dist)) + geom_density(alpha = 0.3) + geom_vline(xintercept = seuil , col = "red")

# ggplot(tab %>% filter(dist == "Grains sans trade off") , aes(sample = ech)) + geom_qq() + geom_qq_line(col = "red")
# 
# ggplot(tab %>% filter(dist == "Grains avec trade off") , aes(sample = ech)) + geom_qq() + geom_qq_line(col = "red")

