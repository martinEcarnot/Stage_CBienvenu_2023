rm(list=ls())

library(ggplot2)
library(tidyverse)
library(shiny)
library(plotly)

# fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


ReRg <- function(NEO , NGO , vg , vinter , vintra , vpos , NGE , nsel){
  NGE * NEO * sqrt(vg+vinter+vpos+vintra) * exp((qnorm(1-nsel/NGO)^2 - qnorm(1-nsel/(NGE*NEO))^2) / 2) / (NGO * sqrt(vg+vinter+vpos+vintra/NGE))
}


isuri <- function(NG_O,NE_O,nsel,NG_E){
  NG_E * NE_O * exp((qnorm(1-nsel/NG_O)^2 - qnorm(1-nsel/(NG_E*NE_O))^2) / 2) / NG_O }

invisuri <- function(NG_O,NE_O,nsel,NG_E){
  1/(NG_E * NE_O * exp((qnorm(1-nsel/NG_O)^2 - qnorm(1-nsel/(NG_E*NE_O))^2) / 2) / NG_O )}






# plot des ie/ig ------------------------------------------------------------


ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(6 , "" ,
           sliderInput("ngs" , "Nombre de grains selectionnes", min = 100 , max = 30000 , step = 1000 , value = 2000),
           sliderInput("nge" , "Nombre de grains par epi", min = 20 , max = 100 , step = 10 , value = 70)),
    
    column(6 , "" ,
           sliderInput("theta","theta" , min = 0 , max = 360 , step = 1 , value = 0),
           sliderInput("phi","phi" , min = -180 , max = 180 , step = 1 , value = 0),
           sliderInput("ngo" , "Nombre de grains observes" , min = 30000 , max = 1000000 , step = 10000 , value = c(70000,500000)),
           sliderInput("neo" , "Nombre d'epis observes" , min = 10 , max = 5000 , step = 10 , value = c(20,1000))
    ),
    fluidRow(
      column(12, "",
             plotOutput("graph"))
    )
  )
)  








server <- function(input, output, session) {
  ngo <- reactive(seq(input$ngo[1] , input$ngo[2] , length = 50))
  neo <- reactive(seq(input$neo[1] , input$neo[2] , length = 50))
  
  t <- reactive(input$theta)
  p <- reactive(input$phi)
  
  sel <- reactive(input$ngs)
  nge <- reactive(input$nge)
  
  zvalue <- reactive(outer(X = ngo() , Y = neo() , FUN = isuri , NG_E = nge() , nsel = sel()) )
  
  couleur <- reactive(apply(zvalue() , MARGIN = c(1,2) , FUN = function(x){
                                                                  if (is.nan(x)==T){return(NA)}
                                                                  if (x<=1){return("#1C86EE")} 
                                                                  if (x>1){return("#32CD32")}
                                                                  }))
 
  
  g <- reactive(persp(x = ngo() ,
                      y = neo() ,
                      z = zvalue() ,
                      theta = t(),
                      phi = p(),
                      xlab = "Nombre de grains observes" ,
                      ylab = "Nombre d'epis observes",
                      zlab = "Repi / Rgrain" , 
                      col = c(couleur())))
  
  
  
  

  output$graph <- renderPlot({g()})
}


shinyApp(ui, server)







# plot des Re/Rg ------------------------------------------------------------


ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(4 , "Autres parametres de selection" ,
           sliderInput("ngs" , "Nombre de grains selectionnes", min = 100 , max = 30000 , step = 1000 , value = 2000),
           sliderInput("nge" , "Nombre de grains par epi", min = 20 , max = 100 , step = 10 , value = 70)),
    
    column(4 , "Axes et parametres graphiques" ,
           sliderInput("theta","theta" , min = 0 , max = 360 , step = 1 , value = 0),
           sliderInput("phi","phi" , min = 0 , max = 360 , step = 1 , value = 0),
           sliderInput("ngo" , "Nombre de grains observes" , min = 40000 , max = 1000000 , step = 20000 , value = c(80000,500000)),
           sliderInput("neo" , "Nombre d'epis observes" , min = 20 , max = 5000 , step = 20 , value = c(20,1000))
          ),
    
    column(4 , "Variances" , 
           sliderInput("vg" , "Vg" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("vintra" , "Vintra" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("vinter" , "Vinter" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("vpos" , "Vpos" , min = 0 , max = 1 , step = 0.1 , value = 1),)
    ),
  
    fluidRow(
      column(12, "",
             plotOutput("graph")))
)





# server <- function(input, output, session) {
#   ngo <- reactive(seq(input$ngo[1] , input$ngo[2] , length = 50))
#   neo <- reactive(seq(input$neo[1] , input$neo[2] , length = 50))
# 
#   t <- reactive(input$theta)
#   p <- reactive(input$phi)
# 
#   sel <- reactive(input$ngs)
#   nge <- reactive(input$nge)
# 
#   vg <- reactive(input$vg)
#   vinter <- reactive(input$vinter)
#   vintra <- reactive(input$vintra)
#   vpos <- reactive(input$vpos)
# 
# 
# 
#   zv <- reactive(outer(X = ngo() , Y = neo() , FUN = isuri , NG_E = nge() , nsel = sel() ))
# 
#   zvalue <- reactive(zv()*sqrt(vg()+vinter()+vpos()+vintra())/sqrt(vg()+vinter()+vpos()+vintra()/nge()))
# 
#   xy <- reactive(expand.grid(ngo(), neo()))
# 
#   g <- {reactive(persp(x = ngo() ,
#                 y = neo() ,
#                 z = zvalue() ,
#                 theta = t(),
#                 phi = p(),
#                 xlab = "Nombre de grains observes" ,
#                 ylab = "Nombre d'epis observes",
#                 zlab = "Repi / Rgrain"))
# 
#     reactive(points(trans3d(xy()[,1], xy()[,2], 1, pmat = g()), col = 2, pch = 10))
#     }
# 
# 
# 
# 
#   output$graph <- renderPlot({g()})
# }


server <- function(input, output, session) {
  ngo <- reactive(seq(input$ngo[1] , input$ngo[2] , length = 50))
  neo <- reactive(seq(input$neo[1] , input$neo[2] , length = 50))
  
  t <- reactive(input$theta)
  p <- reactive(input$phi)
  
  sel <- reactive(input$ngs)
  nge <- reactive(input$nge)
  
  vg <- reactive(input$vg)
  vinter <- reactive(input$vinter)
  vintra <- reactive(input$vintra)
  vpos <- reactive(input$vpos)
  
  
  
  zv <- reactive(outer(X = ngo() , Y = neo() , FUN = isuri , NG_E = nge() , nsel = sel() ))
  
  zvalue <- reactive(zv()*sqrt(vg()+vinter()+vpos()+vintra())/sqrt(vg()+vinter()+vpos()+vintra()/nge()))
  
  plan <- reactive(expand.grid(ngo(), neo()))
  

  
  
  
  
  output$graph <- renderPlot({
    
    graph <- reactive(persp(x = neo() , y = ngo() , z = zvalue() , phi = p() , theta = t()))
    
    for (ix in neo()){
      lines(trans3d(x = ix , y = c(min(ngo()),max(ngo())) , z = 1 , pmat = graph()), col = "red")
      }
    
    for (igrec in ngo()){
      lines(trans3d(x = c(min(neo()),max(neo())) , y = igrec , z = 1 , pmat = graph()), col = "red")
      }
    
      })
}


shinyApp(ui, server)




class(RR(NG_O = 80000 , NE_O = 50 , NG_E = 50 , nsel = 2000 , vinter = 1 , vintra = 1 , vpos = 1 , vg = 1))

ngo <- seq(80000 , 500000 , 20000)

neo <- seq(50 , 1000 , 50)

zvalue <- sapply(ngo , FUN = function(go){sapply(neo , FUN = RR , NG_O = go , NG_E = 70 , nsel = 2000 , vinter = 1 , vintra = 1 , vpos = 1 , vg = 1)}
            )
plan <- expand.grid(neo,ngo)


{

  
  #points(trans3d(plan[,1], plan[,2], 1, pmat = graph), col = "red", pch = 19 , type = "c")
}

# # analytique et tirage donne la même chose ? ------------------------------
# 
# an <- c()
# tir <- c()
# 
# for (vg in seq(0.8 , 1 , 0.1)){
#   for (vinter in seq(0.8 , 1 , 0.1)){
#     for (vintra in seq(0.8 , 1 , 0.1)){
#       for (vpos in seq(0.8 , 1 , 0.1)){
#         for (r in seq(0.8 , 1 , 0.1)){
#           for(NG_E in 68 : 70){
#             for (NE_P in 1:3){
#               for (NE_O in c(4000 , 4500 , 5000)){
#                 
#                 tab <- simulation(mu = mu , 
#                                   vg = vg , 
#                                   vinter = vinter , 
#                                   vintra = vintra , 
#                                   vpos = vpos , 
#                                   r = r , 
#                                   NG_E = NG_E , 
#                                   NE_P = NE_P , 
#                                   NE_O = NE_O)
#                 
#                 s_epi <- sqrt(vg + vinter + vpos)
#                 s_intra <- sqrt(vintra)
#                 vp <- vintra + vg + vinter + vpos
#                 s_p <- sqrt(vp)
#                 p <- r/(NG_E * NE_P)
#                 
#                 # Grains
#                 
#                 grains <- rnorm(n = 1000000 , mean = mu , sd = s_p)
#                 ag <- qnorm(p = 1-p , mean = mu , sd = s_p)
#                 sel <- grains[which(grains>ag)]
#                 Sg2 <- mean(sel) - mu
#                 
#                 # Epis
#                 
#                 epis <- rnorm(n = round(1000000/NG_E) , mean = mu , sd = s_epi)
#                 ae <- qnorm(p = 1-p , mean = mu , sd = s_epi)
#                 epsel <- epis[which(epis > ae)]
#                 grsel <- c()
#                 for (e in epsel){
#                   grsel <- c(grsel , rnorm(n = NG_E , mean = e , sd = s_intra))
#                 }
#                 
#                 Se2 <- mean(grsel) - mu
#                 
#               
#               an <- c(an , tab[1,"Sg"] , tab[1,"Se"])
#               tir <- c(tir , Sg2 , Se2)
#               
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }
# 
# 
# tab <- data.frame(analytique = an , tirage = tir)
# 
# ggplot(data = tab , aes(x = analytique , y = tirage)) + geom_point() + geom_abline(slope = 1 , intercept = 0 , col = "red")
# 
# cor(tir,an) #0.9968019

x <- seq(-10, 10, length= 30)
y <- x
f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "Sinc( r )"
) -> res
round(res, 3)

# (2) Add to existing persp plot - using trans3d() :

xE <- c(-10,10); xy <- expand.grid(xE, xE)
points(trans3d(xy[,1], xy[,2], 6, pmat = res), col = 2, pch = 16)
lines (trans3d(x, y = 10, z = 6 + sin(x), pmat = res), col = 3)

phi <- seq(0, 2*pi, len = 201)
r1 <- 7.725 # radius of 2nd maximum
xr <- r1 * cos(phi)
yr <- r1 * sin(phi)
lines(trans3d(xr,yr, f(xr,yr), res), col = "pink", lwd = 2)
## (no hidden lines)



# (4) Surface colours corresponding to z-values

par(bg = "white")
x <- seq(-1.95, 1.95, length = 30)
y <- seq(-1.95, 1.95, length = 35)
z <- outer(x, y, function(a, b) a*b^2)
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
persp(x, y, z, col = color[facetcol], phi = 30, theta = -30)

par(op)






ngo <- seq(30000,500000,length=50)
neo <- seq(10 , 1000 , length=50)

z <- outer(X = ngo , Y = neo , FUN = isuri , NG_E = 70 , nsel = 2000)
z <- z*sqrt(3+1/70)/sqrt(4)





couleur <- apply(z , MARGIN = c(1,2) , FUN = function(x){
  if (is.nan(x)==T){return(NA)}
  if (x<=1){return("#1C86EE")} 
  if (x>1){return("#32CD32")}
})


a <- persp(x = ngo ,
      y = neo ,
      z = z ,
      xlab = "Nombre de grains observes" ,
      ylab = "Nombre d'epis observes",
      zlab = "Rgrain / Repi",
      col = c(couleur))




xy <- expand.grid(ngo, neo)
points(trans3d(xy[,1], xy[,2], 1, pmat = g), col = 2, pch = 16)










library(plotly)

# en fonction de nge et neo
ngo <- seq(3000,200000,length=100)
neo <- seq(100 , 2000 , length=100)

z <- outer(X = ngo , Y = neo , FUN = invisuri , NG_E = 70 , nsel = 3000)
zvalue <- t(z*sqrt(5+5/70)/sqrt(10))

axy <- list(
  title = "Nombre d'?pis observ?s")

axx <- list(
  title = "Nombre grains observ?s")

axz <- list(
  title = "Rgrain/Repi")

fig <- plot_ly(colors = c('green', 'blue') , showscale = FALSE , x = ~ngo , y = ~neo , z = ~zvalue) %>%
  add_surface() %>% add_surface(z = matrix(1 , ncol = 2 , nrow = 2) , y = c(0,2000) , x = c(0,200000), opacity = 0.5) %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
fig





RR(NE_O = 100,
   NG_O = 500000,
   nsel = 2000,
   NG_E = 70,
   vg = 1,
   vinter = 1,
   vintra = 1,
   vpos = 1)

# en fonction de i/i et h2/h2
rrih <- function(i,h){
  i*sqrt(h)
}

i <- seq(0.1,1,length=50)
h <- seq(1 , 1.5 , length=50)

zvalue <- outer(X = i , Y = h , FUN = rrih)

axy <- list(
  title = "Hepi/Hgrain")

axx <- list(
  title = "iepi/igrain")

axz <- list(
  title = "Repi/Rgrain")

fig <- plot_ly(colors = c('green', 'blue') , showscale = FALSE , x = ~i , y = ~h , z = ~zvalue) %>%
  add_surface() %>% 
  add_surface(z = matrix(1 , ncol = 2 , nrow = 2) , y = c(1,1.5) , x = c(0,1), opacity = 0.5) %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

fig







# brute force pour trouver l'espace paramétrique --------------------------

setwd("~/Stage/Analyses")

load("../donnees/champ")

mean(champ$nb_grain , na.rm = T)
aggregate(champ$nb_grain~champ$selection , FUN = mean , na.rm = T)

load("../donnees/bac")
aggregate(bac$nb_grain ~ bac$BAC , FUN = mean , na.rm = T)


# On fixe NGE à 40

vg <- 1.303^2
vinter <- 1.423^2
vpos <- 1.229^2
vintra <- 2.677^2



# on y va par tranche d'ordre de grandeur
tot <- 88740
param_v_fixe <- data.frame()
i <- 1
for (neo in c(seq(10,100,10),seq(200,1000,100),seq(1250,3000,250),4000,5000)){
  for (ngo in c(seq(10000,100000,10000) , seq(150000,1000000,50000) , seq(1500000,4000000,500000))){
    for (nge in seq(20,100,10)){
      for (nsel in c(200 , 300 , 400 , 1000 , 2000 , 5000 , 10000 , 25000 , 50000 , 100000)){
        
        if (nsel < neo*nge & nsel < ngo){
        param_v_fixe[i,"NEO"] <- neo
        param_v_fixe[i,"NGO"] <- ngo
        param_v_fixe[i,"NGE"] <- nge
        param_v_fixe[i,"nsel"] <- nsel
        param_v_fixe[i,"vg"] <- vg
        param_v_fixe[i,"vintra"] <- vintra
        param_v_fixe[i,"vinter"] <- vinter
        param_v_fixe[i,"vpos"] <- vpos
        
        i <- i+1
        print(round(i*100/tot , 2))
        }
        
      }
    }
  }
}



param_v_fixe$ReRg <- ReRg(NEO = param_v_fixe$NEO , 
                          NGO = param_v_fixe$NGO ,
                          NGE = param_v_fixe$NGE , 
                          nsel = param_v_fixe$nsel , 
                          vg = param_v_fixe$vg , 
                          vinter = param_v_fixe$vinter ,
                          vintra = param_v_fixe$vintra , 
                          vpos = param_v_fixe$vpos)


table(is.na(param_v_fixe$ReRg))

#save(param_v_fixe , file = "../donnees/param_v_fixe")






rm(list = ls())

load("../donnees/param_v_fixe")

library(scales)


donplot3 <- param_v_fixe %>% mutate_at(.vars = c("NEO","nsel","NGE") , .funs = as.factor) %>% mutate(RgRe = 1/ReRg) %>% filter(nsel %in% c(10000,25000,50000,100000)) %>% filter(NGE %in% c(40,60,80)) %>% filter(NEO %in% c(100 , 500 , 1000 , 3000 , 4000 , 5000)) %>% mutate(nsel2 = paste("nsel =",nsel) , NGE2 = paste("NGE =",NGE))

donplot3$nsel3 <- factor(donplot3$nsel2 , levels = c("nsel = 400","nsel = 10000","nsel = 25000","nsel = 50000","nsel = 1e+05"))

ggplot(donplot3 , aes(x = NGO , y = RgRe , col = NEO)) + geom_line(linewidth = 1) + geom_hline(yintercept = 1 , col = "black" , linewidth = 1) + facet_grid(NGE2~nsel3) + theme(axis.text.x = element_text(angle = 45 , hjust = 1) , panel.border = element_rect(colour = "black" , fill=NA)) + labs(y = "R_grain / R_épi" , title = "Valeurs de R_grain / R_épi pour différents jeux de paramètres" , x = "Nombre de grains observés") + ylim(0,3)



ggplot(donplot3 , aes(x = NGO , y = RgRe , col = NEO)) + geom_line(linewidth = 1) + geom_hline(yintercept = 1 , col = "black" , linewidth = 1) + facet_grid(NGE2~nsel3) + theme(axis.text.x = element_text(angle = 45 , hjust = 1) , panel.border = element_rect(colour = "black" , fill=NA) , panel.background = element_blank()) + labs(y = "R_grain / R_epi" , x = "Nombre de grains observ�s") + ylim(0,3) + scale_x_continuous(trans = 'log10')

#+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x) , labels = trans_format("log10", math_format(10^.x)))







# autre situation 

param_v_fixe <- data.frame()
i <- 1
for (neo in seq(2000000,5000000,1000000)){
  for (ngo in seq(100000000 , 3000000000 , 100000000)){
    for (nge in seq(40,90,10)){
      for (nsel in 50000000){
        
        if (nsel < neo*nge & nsel < ngo){
          param_v_fixe[i,"NEO"] <- neo
          param_v_fixe[i,"NGO"] <- ngo
          param_v_fixe[i,"NGE"] <- nge
          param_v_fixe[i,"nsel"] <- nsel
          param_v_fixe[i,"vg"] <- vg
          param_v_fixe[i,"vintra"] <- vintra
          param_v_fixe[i,"vinter"] <- vinter
          param_v_fixe[i,"vpos"] <- vpos
          
          i <- i+1
        }
        
      }
    }
  }
}

param_v_fixe$ReRg <- ReRg(NEO = param_v_fixe$NEO , 
                          NGO = param_v_fixe$NGO ,
                          NGE = param_v_fixe$NGE , 
                          nsel = param_v_fixe$nsel , 
                          vg = param_v_fixe$vg , 
                          vinter = param_v_fixe$vinter ,
                          vintra = param_v_fixe$vintra , 
                          vpos = param_v_fixe$vpos)


param_v_fixe$RgRe <- 1/param_v_fixe$ReRg
param_v_fixe$NEO <- as.factor(param_v_fixe$NEO)

ggplot(param_v_fixe , aes(x = NGO , y = RgRe , col = NEO)) + geom_line(linewidth = 1) + geom_hline(yintercept = 1 , col = "black" , linewidth = 1) + facet_grid(NGE~nsel) + theme(axis.text.x = element_text(angle = 45 , hjust = 1) , panel.border = element_rect(colour = "black" , fill=NA)) + labs(y = "R_grain / R_épi" , title = "Valeurs de R_grain / R_épi pour différents jeux de paramètres" , x = "Nombre de grains observés")









# test d'un graph


param_v_fixe <- data.frame()
i <- 1
for (neo in seq(2000,5000,1)){
  for (ngo in seq(200000,300000,5000)){
    for (nge in 60){
      for (nsel in 100000){
        
        if (nsel < neo*nge & nsel < ngo){
          param_v_fixe[i,"NEO"] <- neo
          param_v_fixe[i,"NGO"] <- ngo
          param_v_fixe[i,"NGE"] <- nge
          param_v_fixe[i,"nsel"] <- nsel
          param_v_fixe[i,"vg"] <- vg
          param_v_fixe[i,"vintra"] <- vintra
          param_v_fixe[i,"vinter"] <- vinter
          param_v_fixe[i,"vpos"] <- vpos
          
          i <- i+1
          print(i)
        }
        
      }
    }
  }
}

param_v_fixe$ReRg <- ReRg(NEO = param_v_fixe$NEO , 
                          NGO = param_v_fixe$NGO ,
                          NGE = param_v_fixe$NGE , 
                          nsel = param_v_fixe$nsel , 
                          vg = param_v_fixe$vg , 
                          vinter = param_v_fixe$vinter ,
                          vintra = param_v_fixe$vintra , 
                          vpos = param_v_fixe$vpos)

param_v_fixe$RgRe <- 1/param_v_fixe$ReRg




ggplot(param_v_fixe , aes(x = NGO , y = RgRe , col = NEO)) + geom_point() + scale_colour_gradientn(colors=rainbow(20)) 
