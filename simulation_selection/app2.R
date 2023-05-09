library(shiny)
library(shinyWidgets)
library(ggplot2)



# Fonctions ---------------------------------------------------------------

i_p <- function(p){
  exp(- qnorm(1-p)^2 / 2) / (p * sqrt(2*pi) )
}

erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

erfc <- function(x) {2 * pnorm(x * sqrt(2), lower = FALSE)}

tronc_epi <- function(x , mu , s_epi , s_intra , seuil){
  -((exp(-((mu-x)^2/(2*(s_epi^2+s_intra^2))))*(s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2+s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2*erf((mu*s_intra^2+s_epi^2*x)/(sqrt(2)*s_epi^2*sqrt(1/s_epi^2+1/s_intra^2)*s_intra^2))-((mu*s_intra^2+s_epi^2*x)*erf(sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/sqrt((mu*s_intra^2+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))+((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)*erf(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2)))/sqrt(2)))/(sqrt((mu*s_intra^2-seuil*(s_epi^2+s_intra^2)+s_epi^2*x)^2/(s_epi^2*s_intra^2*(s_epi^2+s_intra^2))))))/(sqrt(2*pi)*s_epi*s_intra*(s_epi^2+s_intra^2)*(-2+erfc((-seuil+mu)/(sqrt(2)*s_epi)))))
}


maxi <- function(a , m , s_e , s_i){
  x <- a
  y <- a + 0.01
  
  if (is.nan(tronc_epi(x = x , mu = m , s_epi = s_e , s_intra = s_i , seuil = a)) == T){
    x <- y
    y <- y+0.01
  }
  
  while (tronc_epi(x = x , mu = m , s_epi = s_e , s_intra = s_i , seuil = a) < tronc_epi(x = y , mu = m , s_epi = s_e , s_intra = s_i , seuil = a)){
    x <- y
    y <- y + 0.01
  }
  
  (x+y)/2
}


# appli -------------------------------------------------------------------


ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(4 , "Parametres de la parcelle" ,

           #numericInput("coef","Coefficient de changement de surface" , value = 0.5 , step = 0.0001),
           sliderInput("coef","Coefficient de changement de surface" , min = 0.01 , max = 1 , step = 0.01 , value = 0.5),
           sliderInput("surface","Surface initiale" , min = 1 , max = 10000 , step = 1 , value = 100),
           sliderInput("dens","Nombre de plantes par m2" , min = 200 , max = 400 , step = 50 , value = 300)
           ),
    
    column(4 , "Paramètres des plantes",
           sliderInput("grainepi","Nombre de grains par épi" , min = 20 , max = 90 , step = 1 , value = 70),
           sliderInput("epiplante","Nombre d'épis par plante" , min = 1 , max = 5 , step = 1 , value = 2)
           ),
    
    column(4 , "Variances",
           sliderInput("gen","V_G" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("inter","V_inter_epi" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("intra","V_intra_epi" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("pos","V_pos" , min = 0 , max = 1 , step = 0.1 , value = 1)
           )
    ),
  
  fluidRow(
    
    column(4 , "Faisabilité",
           textOutput("p"),
           textOutput("epi_s"),
           textOutput("epi_o")
    ),
    column(4 , "Sélection",
           tableOutput("table")
           )
  ),
  
  fluidRow(
    column(5, "",
           plotOutput("distrib")),
    
    column(2, ""),
    
    column(5, "",
           plotOutput("trunc"))
  )
)

server <- function(input, output, session) {
  
  vg <- reactive(input$gen)
  vpos <- reactive(input$pos)
  vinter <- reactive(input$inter)
  vintra <- reactive(input$intra)
  
  vp <- reactive(vg() + vpos() + vinter() + vintra())
  sp <- reactive(sqrt(vp()))
  
  sinter <- reactive(sqrt(vinter()))
  sintra <- reactive(sqrt(vintra()))
  spos <- reactive(sqrt(vpos()))
  
  sepi <- reactive(sqrt(vg() + vpos() + vinter()))
  
  red  <- reactive(input$coef)
  NG_E <- reactive(input$grainepi)
  NE_P <- reactive(input$epiplante)
  d <- reactive(input$dens)
  surf <- reactive(input$surface)
  
  mu <- 0

  # Calcul de la proportion de grain/epis gardés
  p <- reactive(red() / (NG_E() * NE_P()))
  
  
  # Modelisation grain a grain 
  ag <- reactive(
    qnorm(p = 1 - p() , mean = mu , sd = sp())    )
  
  ig <- reactive(
    i_p(p())    )
  
  H2g <- reactive(
    vg() / vp()  )
  
  Sg <- reactive(
    ig() * sp()  )
  
  Rg <- reactive(
    Sg() * H2g()
  )
  
  
  
  # Modelisation epi par epi
  
  ae <- reactive(
    qnorm(p = 1-p() , mean = mu , sd = sepi())  )
  
  Se <- reactive(
    maxi(a = ae() , m = mu , s_e = sepi() , s_i = sintra())  )
  
  ie <- reactive(
    Se() / sp()  )
  
  H2e <- reactive(
    vg() / (vg() + vpos() + vinter() + vintra()/NG_E())  )
  
  Re <- reactive(
    H2e() * Se()  )
  
  NE_S <- reactive(           # Nombre d'épis sélectionnés
    p() * NE_P() * d() * surf()  )
  
  NE_O <- reactive(           # nombre d'épis dans la parcelle
    NE_P() * d() * surf()  )
      
  
  
  
  # Outputs
  
  ## tableau des résultats
  res <- reactive(
    data.frame(row.names = c("Sélection sur grain" , "Sélection sur épi"),
               H2 = c(round(H2g() , 2) , round(H2e() , 2)),
               i = c(round(ig() , 2) , round(ie() , 2)),
               R = c(round(Rg() , 2) , round(Re() , 2)),
               S = c(round(Sg() , 2) , round(Se() , 2)))
    )
  
  
  # Graph des distributions des grains et des moyennes des épis
  x <- reactive(
    seq(mu - 3.5 * sp() , mu + 3.5 * sp() , 0.1)    )

  don <- reactive(
    data.frame(abs = rep(x(),2) ,
               ord =c(dnorm(x = x() , mean = mu , sd = sp()) , dnorm(x = x() , mean = mu , sd = sepi())) ,
               Distribution = c(rep("Grains" , length(x())) , rep("Moyenne des épis" , length(x()))))
    )
  
  p1 <- reactive(
    ggplot(data = don() , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = mu) + labs(x = "Phenotype du grain" , y = "Densité de probabilité" , title = "Distributions des grains et des épis") + geom_vline(xintercept = ag() , col = "#F8766D") + geom_vline(xintercept = ae() , col = "#619CFF") + scale_colour_manual(values = c("#F8766D", "#619CFF")) + theme(plot.background = element_rect(colour = "black"))
  )


  
  ## Graph des distributions tronquées
  
  abs_g <- reactive(
    seq(ag() , mu + 4 * sp() , 0.1))
  
  abs_e <- reactive(
    seq(Se() + mu - 2 * sintra() , Se() + mu + 2 * sintra() , 0.1))
  
  don2 <- reactive(
    data.frame(abs = c(abs_g() , abs_e()) ,
               ord = c(dnorm(x = abs_g() , mean = mu , sd = sp()) / (1-pnorm(q = ag() , mean = mu , sd = sp())) , tronc_epi(x = abs_e() , mu = mu , s_epi = sepi() , s_intra = sintra() , seuil = ae())) ,
               Distribution = c(rep("Sélection sur grain" , length(abs_g())) , rep("Sélection sur épi" , length(abs_e()))))
    )
  
  
  p2 <- reactive(
    ggplot(data = don2() , aes(x = abs , y = ord , col = Distribution)) + geom_line() + geom_vline(xintercept = ag() , col = "#F8766D") + geom_vline(xintercept = ae() , col = "#619CFF") + scale_colour_manual(values = c("#619CFF", "#F8766D")) + labs(x = "Phenotype du grain" , y = "Densité de probabilité" , title = "Distributions des grains sélectionnés") + theme(plot.background = element_rect(colour = "black"))
    )
  
  
  # Lien avec ui
  output$table <- renderTable({res()} , rownames = TRUE)
  
  output$p <- renderText({paste0("P% = " , round(p()*100 , 6))})
  output$epi_s <- renderText({paste0("Nombre d'épis à sélectionner = " , round(NE_S()))})
  output$epi_o <- renderText({paste0("Nombre d'épis dans la parcelle = " , round(NE_O()))})
  
  output$distrib <- renderPlot({p1()})
  output$trunc <- renderPlot({p2()})
  
}

shinyApp(ui, server)