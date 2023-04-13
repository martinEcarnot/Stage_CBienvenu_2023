rm(list=ls())


# Paramètres

sigma_g <- 1
sigma_e <- 2
sigma_inter <- 0.5
sigma_intra <- 0.3
red  <- 0.9
nb_grain_epi <- 70
nb_epi_plante <- 2
dens <- 300

sigma_p <- sigma_g + sigma_e + sigma_inter + sigma_intra

# Calcul de pourcentage de grains gardés
p <- red/(nb_grain_epi * nb_epi_plante)

# Modelisation grain a grain ----------------------------------------------

# calcul du i correspondant
grains <- rnorm(n = nb_grain_epi * nb_epi_plante * dens , 
                mean = 0 , 
                sd = sigma_p)

seuil_grains <- qnorm(p = 1 - p , 
               mean = 0 , 
               sd = sigma_p)

grains_sel_grains <- grains[which(grains > seuil_grains)]

S_grains <- mean(grains_sel_grains)

i_grains <- S_grains/sigma_p

h2_grains <- sigma_g^2/sigma_p^2

R_grains <- h2_grains*i_grains*sigma_p



# Simulation epi par epi --------------------------------------------------

epis <- rnorm(n = nb_epi_plante * dens , 
              mean = 0 , 
              sd = sigma_inter)

seuil_epis <- qnorm(p = 1 - p , 
               mean = 0 , 
               sd = sigma_inter)

epis_sel <- epis[which(epis > seuil_epis)]

grains_sel_epis <- as.vector(
  sapply(epis_sel , FUN = rnorm , n = nb_grain_epi , sd = sigma_intra))

S_epis <- mean(grains_sel_epis)

i_epis <- S_epis/sigma_p

h2_epis <- sigma_g^2/(sigma_g + sigma_e + sigma_inter + sigma_intra/nb_grain_epi)

R_epis <- h2_epis*i_epis*sigma_p
