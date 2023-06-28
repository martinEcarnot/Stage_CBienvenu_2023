library(shiny)

ui <- fluidPage(
  
  # parametres
  fluidRow( # faire une ligne de parametres
    
    column(6 , "Parametres de la parcelle" ,
           numericInput("coef","Coefficient de changement de surface" , value = 0.5 , step = 0.0001),
           
          # sliderInput("coef","Coefficient de changement de surface" , min = 0 , max = 1 , step = 0.001 , value = 1),
           sliderInput("grainepi","Nombre de grains par épi" , min = 20 , max = 90 , step = 1 , value = 70),
           sliderInput("epiplante","Nombre d'épis par plante" , min = 1 , max = 5 , step = 1 , value = 2),
           sliderInput("dens","Nombre de plantes par m2" , min = 200 , max = 400 , step = 50 , value = 300)
           ),
    
    column(6 , "Variances",
           sliderInput("gen","V_G" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("inter","V_inter_epi" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("intra","V_intra_epi" , min = 0 , max = 1 , step = 0.1 , value = 1),
           sliderInput("env","V_pos" , min = 0 , max = 1 , step = 0.1 , value = 1)
           )
    ),
  
  fluidRow(
    
    column(6 , "Sélection grain à grain",
           textOutput("Hg"),
           textOutput("ig"),
           textOutput("Sg"),
           textOutput("Rg")
           ),
    
    column(6 , "Sélection epi par epi",
           textOutput("He"),
           textOutput("ie"),
           textOutput("Se"),
           textOutput("Re")
    ),
    
    fluidRow(
      column(12 , "P%" , 
             textOutput("p"))
    )
  )
)

server <- function(input, output, session) {
  
  vg <- reactive(input$gen)
  ve <- reactive(input$env)
  vinter <- reactive(input$inter)
  vintra <- reactive(input$intra)
  red  <- reactive(input$coef)
  nb_grain_epi <- reactive(input$grainepi)
  nb_epi_plante <- reactive(input$epiplante)
  dens <- reactive(input$dens)
  
  vp <- reactive(vg() + ve() + vinter() + vintra())
  sigma_p <- reactive(sqrt(vp()))
  
  # Calcul de la proportion de grain/epis gardés
  p <- reactive(red() / (nb_grain_epi() * nb_epi_plante()))
  
  
  # Modelisation grain a grain 
  grains <- reactive(
    rnorm(n = nb_grain_epi() * nb_epi_plante() * dens() , mean = 0 , sd = sigma_p()))
  
  seuil_grains <- reactive(
    qnorm(p = 1 - p() , mean = 0 , sd = sigma_p()))
  
  grains_sel_grains <- reactive(
    grains()[which(grains() > seuil_grains())])
  
  S_grains <- reactive(
    mean(grains_sel_grains()))
  
  i_grains <- reactive(
    S_grains()/sigma_p())
  
  h2_grains <- reactive(
    vg()/vp())
  
  R_grains <- reactive(
    h2_grains() * i_grains() * sigma_p())
  
  
  # Modelisation epi par epi

  sigma_intra <- reactive(
    sqrt(vintra()))
  
  sigma_epi <- reactive(
    sqrt(vg() + ve() + vinter() )
  )
  
  epis <- reactive(
    rnorm(n = nb_epi_plante() * dens() , mean = 0 , sd = sigma_epi() ))
  
  seuil_epis <- reactive(
    qnorm(p = 1 - p() , mean = 0 , sd = sigma_epi() ))
  
  epis_sel <- reactive(
    epis()[which(epis() > seuil_epis())])
  
  grains_sel_epis <- reactive(
    as.vector(sapply(epis_sel() , FUN = rnorm , n = nb_grain_epi() , sd = sigma_intra())) )
  
  S_epis <- reactive(
    mean(grains_sel_epis()))
  
  i_epis <- reactive(
    S_epis() / sigma_p())
  
  h2_epis <- reactive(
    vg()/(vg() + ve() + vinter() + vintra()/nb_grain_epi()))
  
  R_epis <- reactive(
    h2_epis() * i_epis() * sigma_p())
  
  output$Hg <- renderText({paste0("H2 = " , round(h2_grains() , 2))})
  output$ig <- renderText({paste0("i = " , round(i_grains() , 2))})
  output$Sg <- renderText({paste0("S = " , round(S_grains() , 2))})
  output$Rg <- renderText({paste0("R = " , round(R_grains() , 2))})
  
  
  output$He <- renderText({paste0("H2 = " , round(h2_epis() , 2))})
  output$ie <- renderText({paste0("i = " , round(i_epis() , 2))})
  output$Se <- renderText({paste0("S = " , round(S_epis() , 2))})
  output$Re <- renderText({paste0("R = " , round(R_epis() , 2))})
  output$p <- renderText({paste0("P% = " , round(p()*100 , 6))})
}

shinyApp(ui, server)