rm(list=ls())


load("../donnees/bac")



don <- bac %>% filter(is.na(nb_epi)==F & nb_epi > 1 & geno != "INCONNU")

hist(don$nb_epi)


length(unique(don$geno))



x9y3 <- don %>% filter(BAC2 == "x09y03") %>% select(geno , nb_epi) %>% na.omit()

x9y4 <- don %>% filter(BAC2 == "x09y04") %>% select(geno , nb_epi) %>% na.omit()

x9y5 <- don %>% filter(BAC2 == "x09y05") %>% select(geno , nb_epi) %>% na.omit()

x10y3 <- don %>% filter(BAC2 == "x10y03") %>% select(geno , nb_epi) %>% na.omit()

x10y4 <- don %>% filter(BAC2 == "x10y04") %>% select(geno , nb_epi) %>% na.omit()

x10y5 <- don %>% filter(BAC2 == "x10y05") %>% select(geno , nb_epi) %>% na.omit()


a <- merge(x9y3,x9y4 , by = "geno")

b <- merge(x9y3,x9y5 , by = "geno")

c <- merge(x9y3,x10y3 , by = "geno")

d <- merge(x9y3,x10y4 , by = "geno")

e <- merge(x9y3,x10y5 , by = "geno")

f <- merge(x9y4,x9y5 , by = "geno")

g <- merge(x9y4,x10y3 , by = "geno")

h <- merge(x9y4,x10y4 , by = "geno")

i <- merge(x9y4,x10y5 , by = "geno")

j <- merge(x9y5,x10y3 , by = "geno")

k <- merge(x9y5,x10y4 , by = "geno")

l <- merge(x9y5,x10y5 , by = "geno")

m <- merge(x10y3,x10y4 , by = "geno")

n <- merge(x10y3,x10y5 , by = "geno")

o <- merge(x10y4,x10y5 , by = "geno")


a2 <- sample(unique(a$geno) , 7)
geno_ech <- a2

b2 <- ifelse(b$geno %in% geno_ech , NA , b$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , b2)

c2 <- ifelse(c$geno %in% geno_ech , NA , c$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , c2)

d2 <- ifelse(d$geno %in% geno_ech , NA , d$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , d2)

e2 <- ifelse(e$geno %in% geno_ech , NA , e$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , e2)

f2 <- ifelse(f$geno %in% geno_ech , NA , f$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , f2)

g2 <- ifelse(g$geno %in% geno_ech , NA , g$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , g2)

h2 <- ifelse(h$geno %in% geno_ech , NA , h$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , h2)

i2 <- ifelse(i$geno %in% geno_ech , NA , i$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , i2)

j2 <- ifelse(j$geno %in% geno_ech , NA , j$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , j2)

k2 <- ifelse(k$geno %in% geno_ech , NA , k$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , k2)

l2 <- ifelse(l$geno %in% geno_ech , NA , l$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , l2)

m2 <- ifelse(m$geno %in% geno_ech , NA , m$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , m2)

n2 <- ifelse(n$geno %in% geno_ech , NA , n$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , n2)

o2 <- ifelse(o$geno %in% geno_ech , NA , o$geno) %>% na.omit() %>% sample(7)
geno_ech <- c(geno_ech , o2)





length(unique(geno_ech))




ech_x9y3 <- c(a2,b2,c2,d2,e2)

ech_x9y4 <- c(a2,f2,g2,h2,i2)

ech_x9y5 <- c(b2,f2,j2,k2,l2)

ech_x10y3 <- c(c2,g2,j2,m2,n2)

ech_x10y4 <- c(d2,h2,k2,m2,o2)

ech_x10y5 <- c(e2,i2,l2,n2,o2)



table(ech_x9y3 %in% ech_x9y4)

table(ech_x9y3 %in% ech_x9y5)

table(ech_x9y3 %in% ech_x10y3)

table(ech_x9y3 %in% ech_x10y4)

table(ech_x9y3 %in% ech_x10y5)

table(ech_x9y4 %in% ech_x9y5)

table(ech_x9y4 %in% ech_x10y3)

table(ech_x9y4 %in% ech_x10y4)

table(ech_x9y4 %in% ech_x10y5)

table(ech_x9y5 %in% ech_x10y3)

table(ech_x9y5 %in% ech_x10y4)

table(ech_x9y5 %in% ech_x10y5)

table(ech_x10y3 %in% ech_x10y4)

table(ech_x10y3 %in% ech_x10y5)

table(ech_x10y4 %in% ech_x10y5)


# youpi !!
