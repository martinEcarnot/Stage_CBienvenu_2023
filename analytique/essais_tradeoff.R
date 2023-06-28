library(ggplot2)
library(tidyverse)

rm(list = ls())

vg <- 1
vintra <- 5
ve <- 2

NGE <- 70

sepi <- sqrt(vg+ve+vintra/NGE)

p <- 0.2

pval <- S_sans <- S_avec <- c()



for (i in 1:100){

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
a <- -7
oo <- 70

# population de grains

grains_sel_to <- c()

for (e in epis_sel){
  grains_sel_to <- c(grains_sel_to , rnorm(mean = e , sd = vintra , n = oo + a*e))
}

pval <- c(pval,t.test(grains_sel,grains_sel_to)$p.value)

S_sans <- c(S_sans , mean(grains_sel))

S_avec <- c(S_avec , mean(grains_sel_to))
}


plot(S_sans , S_avec)
abline(a = 0 , b = 1)
lm(S_sans ~ S_avec)

hist(pval)

tab <- data.frame(ech = c(epis , grains_sel , grains_sel_to) , 
                  dist = c(rep("moyenne épis" , length(epis)) , rep("Grains sans trade off" , length(grains_sel)) , rep("Grains avec trade off" , length(grains_sel_to))))

moy <- as.data.frame(tab %>% group_by(dist) %>% summarise(moyenne = mean(ech)))

ggplot(tab , aes(x = ech , col = dist , fill = dist)) + geom_density(alpha = 0.3) + geom_vline(xintercept = moy[1,2] , col = "red") + geom_vline(xintercept = moy[2,2] , col = "green") 








# tests formule avec trade off --------------------------------------------

rm(list=ls())

# trade off : NGE = b - a*trait

### simulation
b <- 70
a <- 5
sp <- 3
NEO <- 100
mu <- 5
# distrib epis avec sigma 
epis <- rnorm(n=NEO , sd = sp , mean=mu)

nb_grains <- 0
for (e in epis){
  nb_grains <- nb_grains + b - a*e
}


f <- function(x){x*(b-a*x)}

integrate(f = f , lower = mu-3*sp , upper = mu+3*sp)

plot(seq(mu-3*sp ,mu+3*sp,0.1) , f(seq(mu-3*sp ,mu+3*sp,0.1)))
