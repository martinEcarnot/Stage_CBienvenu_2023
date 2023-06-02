rm(list=ls())

library(ggplot2)
library(shiny)
library(plot3D)
library(tidyverse)

# fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}


RR <- function(vg , vinter , vintra , vpos , NG_E , NE_O , NG_O ,nsel){
  NG_E * NE_O * sqrt(vg+vinter+vpos+vintra) * exp((qnorm(1-nsel/NG_O)^2 - qnorm(1-nsel/(NG_E*NE_O))^2) / 2) / (NG_O * sqrt(vg+vinter+vpos+vintra/NG_E))
}


isuri <- function(nsel,NG_O,NE_O,NG_E){
  NG_E * NE_O * exp((qnorm(1-nsel/NG_O)^2 - qnorm(1-nsel/(NG_E*NE_O))^2) / 2) / NG_O }



simulation <- function(mu , vg , vinter , vintra , vpos , r , NG_E , NE_P , NE_O){
  
  # calculs utiles
  s_epi <- sqrt(vg + vinter + vpos + vintra/NG_E)
  vp <- vintra + vg + vinter + vpos
  s_p <- sqrt(vp)
  p <- r/(NG_E * NE_P)
  
  # Selection grain a grain 
  ig <- i_p(p)
  H2g <- vg / vp
  Sg <- ig * s_p
  Rg <- H2g*Sg
  
  # Selection epi par epi 
  Se <- ig*s_epi
  ie <- Se/s_p
  H2e <- vg / (vg + vinter + vpos + vintra/NG_E)
  Re <- H2e*Se
  NE_S <- NE_O * p
  
  
  # rapport des R
  
  RR <- s_p / s_epi 
  
  # resultats
  data.frame(mu=mu, vg=vg, vinter=vinter , vintra=vintra , vpos=vpos , r=r , NG_E=NG_E , NE_P=NE_P , NE_O=NE_O , p=p ,
             H2g=H2g , ig=ig , Sg=Sg , Rg=Rg ,
             H2e=H2e , ie=ie , Se=Se , Re=Re ,
             RR=RR)
  

}





# plot des i/i ------------------------------------------------------------




ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(6 , "" ,
           sliderInput("Nombre de grains sélectionnés","ngs" , min = 1000 , max = 1000000 , step = 1000 , value = 500000),
           sliderInput("Nombre de grains par épi","nge" , min = 20 , max = 100 , step = 10 , value = 70)),
    
    column(6 , "" ,
           sliderInput("theta","theta" , min = 0 , max = 360 , step = 1 , value = 0),
           sliderInput("phi","phi" , min = 0 , max = 360 , step = 1 , value = 0)
    ),
    fluidRow(
      column(12, "",
             plotOutput("graph"))
    )
  )
)  



ng_o <- seq(70000,1000000,1000)
ne_o <- seq(100,1000,50)
ngo <- c()
neo <- c()
for (a in ng_o){
  for (b in ne_o){
    ngo <- c(ngo,a)
    neo <- c(neo,b)
  }
}

rm(ng_o,ne_o,a,b)

isuri(nsel = 1000000 , NG_O = ngo , NE_O = neo , NG_E = 20)



server <- function(input, output, session) {
  t <- reactive(input$theta)
  p <- reactive(input$phi)
  
  nsel <- reactive(input$ngs)
  nge <- reactive(input$nge)
  
  
  for (a in seq(nsel()+1000,1000000,1000)){
    for (b in seq(round(nsel()/NG_E())+50,5000,50)){
      ngo <- c(ngo,a)
      neo <- c(neo,b)
      ii <- reactive(c(ii(),isuri(nsel=nsel(),NG_O=a,NE_O=b,NG_E=nge())))
    }
  }
  
  g <- reactive(scatter3D(x = neo , y = ngo , z = ii , theta = t() , phi = p())
  )
  
  output$graph <- renderPlot({g()})
}


shinyApp(ui, server)






# plot des R/R en fct de i/i et vintra/vp ---------------------------------

rm(list=ls())

RRii <- function(vg , vinter , vintra , vpos , ii , NG_E){
  ii * sqrt(vg+vinter+vpos+vintra)  / sqrt(vg+vinter+vpos+vintra/NG_E)
}

i <- c()
v <- c()
r <- c()

for (ii in seq(0.01,2,0.01)){
  for (vintra in seq(0.1,4,0.1)){
    v <- c(v , vintra/(vintra+3))
    i <- c(i , ii)
    r <- c(r , RRii(vg=1,vinter=1,vintra=vintra,vpos=1,NG_E=70,ii=ii))
  }
}

ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(12 , "" ,
           sliderInput("theta","theta" , min = 0 , max = 360 , step = 1 , value = 0),
           sliderInput("phi","phi" , min = 0 , max = 360 , step = 1 , value = 0)
    ),
    fluidRow(
      column(12, "",
             plotOutput("graph"))
    )
  )
)  


server <- function(input, output, session) {
  t <- reactive(input$theta)
  p <- reactive(input$phi)
  
  g <- reactive(scatter3D(x = i , y = v , z = r , theta = t() , phi = p())
  )
  
  output$graph <- renderPlot({g()})
}


shinyApp(ui, server)



# # analytique et tirage donne la mÃªme chose ? ------------------------------
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

