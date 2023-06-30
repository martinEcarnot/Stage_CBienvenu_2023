rm(list=ls())

setwd("~/Stage/Analyses")

load("../donnees/bac")

og <- bac

bac <- og %>% filter(semis == "06/01" & is.na(preco)==F) %>% mutate(preco_voisin = NA)

for (i in 1:nrow(bac)){
  
  # coordonnees de la plante focale
  xf <- bac[i,"X"]
  yf <- bac[i,"Y"]
  bacf <- bac[i,"BAC"]
  
  # precocite des plantes voisines
  p1 <- bac[which(bac$X == xf+1 & bac$Y == yf & bac$BAC == bacf),"preco"]
  p2 <- bac[which(bac$X == xf-1 & bac$Y == yf & bac$BAC == bacf),"preco"]
  p3 <- bac[which(bac$X == xf & bac$Y == yf+1 & bac$BAC == bacf),"preco"]
  p4 <- bac[which(bac$X == xf & bac$Y == yf-1 & bac$BAC == bacf),"preco"]
  p5 <- bac[which(bac$X == xf+1 & bac$Y == yf+1 & bac$BAC == bacf),"preco"]
  p6 <- bac[which(bac$X == xf-1 & bac$Y == yf-1 & bac$BAC == bacf),"preco"]
  p7 <- bac[which(bac$X == xf+1 & bac$Y == yf-1 & bac$BAC == bacf),"preco"]
  p8 <- bac[which(bac$X == xf-1 & bac$Y == yf+1 & bac$BAC == bacf),"preco"]
  
  # on garde que ceux qui sont pas vides
  la <- c()
  
  for (p in c(p1,p2,p3,p4,p5,p6,p7,p8)){
    if (is_empty(p)==F){la <- c(la,p)}
  }
  
  # On fait la moyenne et on la met dans le tableau
  bac[i,"preco_voisin"] <- mean(la)
  
}


rm(p1,p2,p3,p4,p5,p6,p7,p8,la,i,p,xf,yf,bacf)

ggplot(data=bac , aes(x = preco , y=preco_voisin)) + geom_point() + geom_smooth(method = "lm")

ggplot(data=bac , aes(x = preco , y=preco_voisin)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~BAC)

ggplot(data=bac , aes(x = preco , y=preco_voisin)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~luz)

ggplot(data=bac , aes(x = preco , y=preco_voisin)) + geom_point() + geom_smooth(method = "lm") + facet_grid(BAC~luz)


cor(bac$preco,bac$preco_voisin)

mod <- lm(preco ~ preco_voisin , data = bac)
par(mfrow = c(2,2))
plot(mod)
summary(mod)

mod <- lm(preco_voisin ~ preco + BAC , data = bac)
par(mfrow = c(2,2))
plot(mod)
summary(mod)

mod <- lm(preco_voisin ~ preco * BAC , data = bac)
par(mfrow = c(2,2))
plot(mod)
summary(mod)

mod <- lm(preco_voisin ~ preco + luz, data = bac)
par(mfrow = c(2,2))
plot(mod)
summary(mod)

mod <- lm(preco_voisin ~ preco + BAC + luz, data = bac)
par(mfrow = c(2,2))
plot(mod)
summary(mod)



ggplot(bac , aes(x = X , y = Y , fill = preco)) + geom_tile() + facet_wrap(~BAC)
ggplot(bac , aes(x = X , y = Y , fill = preco_voisin)) + geom_tile() + facet_wrap(~BAC)
ggplot(bac , aes(x = X , y = Y , fill = preco - preco_voisin)) + geom_tile() + facet_wrap(~BAC)


ggplot(data=bac , aes(x = preco , y=preco_voisin)) + geom_point()
