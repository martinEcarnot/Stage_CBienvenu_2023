library(ggplot2)
library(tidyverse)


i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

i_p(1)
i_p(0.5)
i_p(0.07)
i_p(0.26)
i_p(0.3)
i_p(0.6)




i <- sapply(seq(0.01,1,0.01) , FUN = i_p)

plot(x = seq(0.01,1,0.01) , y = i)




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



# calcul intégrale par méthode numérique ----------------------------------


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
se <- 1.5
si <- 0.5
a <- 10

abs <- seq(mu_epi - 4 * se , mu_epi + 5 * se , 0.01)

dist_epi <- sapply(abs , FUN = dnorm , mean = mu_epi , sd = se)

dist_grain <- sapply(abs , FUN = f , 
                     seuil = a,
                     mepi = mu_epi,
                     epi = se,
                     intra = si)

dist_epi <- as.data.frame(dist_epi) %>% rename("Y" = "dist_epi") %>% mutate(distribution = "EPI" , X = abs)

don <- as.data.frame(dist_grain) %>% rename("Y" = "dist_grain") %>% mutate(distribution = "GRAINS" , X = abs) %>% rbind(dist_epi)

ggplot(data = don , aes(x = X , y = Y , col = distribution )) + geom_line() + geom_vline(xintercept = a)

}

 qqnorm(dist_grain)
 