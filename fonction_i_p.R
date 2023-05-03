rm(list = ls())

library(ggplot2)
library(tidyverse)


i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


i <- sapply(seq(0.01,1,0.01) , FUN = i_p)

plot(x = seq(0.01,1,0.01) , y = i)
rm(i)


# Error function

erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

erfc <- function(x) {2 * pnorm(x * sqrt(2), lower = FALSE)}


tronc_epi <- function(x , mue , se , si , a){
  -((exp(-((mue-x)^2/(2*(se^2+si^2))))*(se^2*sqrt(1/se^2+1/si^2)*si^2+se^2*sqrt(1/se^2+1/si^2)*si^2*erf((mue*si^2+se^2*x)/(sqrt(2)*se^2*sqrt(1/se^2+1/si^2)*si^2))-((mue*si^2+se^2*x)*erf(sqrt((mue*si^2+se^2*x)^2/(se^2*si^2*(se^2+si^2)))/sqrt(2)))/sqrt((mue*si^2+se^2*x)^2/(se^2*si^2*(se^2+si^2)))+((mue*si^2-a*(se^2+si^2)+se^2*x)*erf(sqrt((mue*si^2-a*(se^2+si^2)+se^2*x)^2/(se^2*si^2*(se^2+si^2)))/sqrt(2)))/(sqrt((mue*si^2-a*(se^2+si^2)+se^2*x)^2/(se^2*si^2*(se^2+si^2))))))/(sqrt(2*pi)*se*si*(se^2+si^2)*(-2+erfc((-a+mue)/(sqrt(2)*se)))))
}



# parametres de la simulation ---------------------------------------------
mu_epi <- 5
s_epi <- 3
s_intra <- 2
s_env <- 7
s_g <- 4
seuil <- 8

abs <- seq(mu_epi-10 , mu_epi + 10 * s_epi , 0.1)

# fonction exacte troncation epi ------------------------------------------
{
dist_epi <- sapply(abs , FUN = dnorm , mean = mu_epi , sd = s_epi)

dist_grain <- sapply(abs , FUN = tronc_epi , 
                     a = seuil,
                     mue = mu_epi,
                     se = s_epi,
                     si = s_intra)

dist_epi <- as.data.frame(dist_epi) %>% rename("Y" = "dist_epi") %>% mutate(distribution = "EPI" , X = abs)

don <- as.data.frame(dist_grain) %>% rename("Y" = "dist_grain") %>% mutate(distribution = "GRAINS" , X = abs) %>% rbind(dist_epi) %>% na.omit()



p <- ggplot(data = don , aes(x = X , y = Y , col = distribution )) + geom_line() + geom_vline(xintercept = seuil)

# est-ce que tron epi est bien une fonction de densité : 

I <- integrate(tronc_epi , lower = -Inf , upper = Inf ,
          a = seuil,
          mue = mu_epi,
          se = s_epi,
          si = s_intra)
# Oui ça s'intègre en 1


# Est-ce que c'est une densité normale ?
# tirage dans la loi :
tab <- subset(don , distribution == "GRAINS" , select = c("X","Y"))

tirage <- c()
for (i in 1:nrow(tab)){
  tirage <- c(tirage , rep(tab[i,1] , round(tab[i,2] * 1000)) )
}

qqnorm(tirage)
print(p)
I

rm(I,p,i,tirage,dist_epi,dist_grain,don,tab)
}



# calcul int?grale par m?thode num?rique ----------------------------------


integ <- function(t,ix,s,mue,sigma_e,sigma_i){
  1/(1-pnorm(q = s , mean = mue , sd = sigma_e)) *   
    dnorm(x = t , mean = mue , sd = sigma_e)  *  dnorm(x = ix , mean = t , sd = sigma_i)
}



f <- function(x,seuil,mepi,epi,intra){
  integrate(f = integ , lower = seuil , upper =  Inf,
            ix = x,
            s = seuil,
            mue = mepi,
            sigma_e = epi,
            sigma_i = intra)$value
}

{
  mu_epi <- 6
  se <- s_epi
  si <- s_intra
  a <- seuil
  
    dist_epi <- sapply(abs , FUN = dnorm , mean = mu_epi , sd = se)
  
  dist_grain <- sapply(abs , FUN = f , 
                       seuil = a,
                       mepi = mu_epi,
                       epi = se,
                       intra = si)
  
  dist_epi <- as.data.frame(dist_epi) %>% rename("Y" = "dist_epi") %>% mutate(distribution = "EPI" , X = abs)
  
  don <- as.data.frame(dist_grain) %>% rename("Y" = "dist_grain") %>% mutate(distribution = "GRAINS" , X = abs) %>% rbind(dist_epi)
  
  print(ggplot(data = don , aes(x = X , y = Y , col = distribution )) + geom_line() + geom_vline(xintercept = a))
  
  rm(a,dist_grain,se,si,dist_epi,don)
}

# numérique et analytique ça donne la même chose



# verif que integrer distrib epi donne bien distrib grains ----------------

s <- sqrt(s_epi^2 + s_g^2 + s_env^2)

s_grain <- sqrt(s_epi^2 + s_g^2 + s_env^2 + s_intra^2)

dist_epi <- dnorm(x = abs , mean = mu_epi , sd = s)

dist_grain <- dnorm(x = abs , mean = mu_epi , sd = s_grain)


integ <- function(t,ix,s,mue,sigma_e,sigma_i){
  dnorm(x = t , mean = mue , sd = sigma_e)  *  dnorm(x = ix , mean = t , sd = sigma_i)
}



f <- function(x,mepi,epi,intra){
  integrate(f = integ , lower = -Inf , upper =  Inf,
            ix = x,
            s = seuil,
            mue = mepi,
            sigma_e = epi,
            sigma_i = intra)$value
}

dist_integ <- sapply(abs , FUN = f ,
                     mepi = mu_epi,
                     epi = s,
                     intra = s_intra)



don <- data.frame(X = rep(abs,3),
                  Y = c(dist_grain,dist_epi,dist_integ),
                  distribution = c(rep("grain",length(abs)) , rep("epi",length(abs)) , rep("integrale",length(abs))))


ggplot(data = don , aes(x = X , y = Y , col = distribution)) + geom_line()



# Putain ça marche !!!!!



# tests -------------------------------------------------------------------

mu <- 5
sigma <- 0.5

a <- mu + 2*sigma

absciss <- seq(mu-4*sigma,mu+4*sigma,0.01)
n <- sapply(absciss , FUN = dnorm , mean = mu , sd = sigma)

p <- sapply(absciss , FUN = pnorm , mean = a , sd = 0.001)


ggplot() + 
  geom_area(aes(x = absciss , y = n) , fill = "blue" , alpha = 0.3) + 
  geom_line(aes(x = absciss , y = p) , col = "red") +
  geom_area(aes(x = absciss , y = n*p) , fill = "green" , alpha = 0.3)


