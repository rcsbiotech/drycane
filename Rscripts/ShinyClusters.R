# TEST PAGE FOR SHINY APP #
# MODIFY FREELY BEFORE SAVING ON TAB 3 #

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Cana-Seca, Tempo e Estresse"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("CLA",
                  "Cluster A: ",
                  min = 1,
                  max = 9,
                  value = 1),
      sliderInput("CLB",
                  "Cluster B: ",
                  min = 1,
                  max = 9,
                  value = 2),
      sliderInput("DaysOfTime",
                  "Tempo (30, 60, 90): ",
                  min = 30,
                  max = 90,
                  value = 30,
                  step = 30),
      selectInput("TagPicker", "Isoforma: ",
                  as.data.frame(AmbosClustersSig$at30, stringsAsFactors=F)$Tag)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    # Input01: Qual data DEG (30, 60, 90)
    # Diferencialmente expresso na data..
    InputDate <- input$DaysOfTime
    
    # Qual tag consultar?
    QueryTag = input$TagPicker
    
    ## Filtrando os valores para cada tag
    # Em 30 dias
    Tag30Tolerant <- subset(Sig30, (Tag == QueryTag & Genotype == "T"))
    Tag30Susceptible <- subset(Sig30, (Tag == QueryTag & Genotype == "S"))
    
    # Em 60 dias
    Tag60Tolerant <- subset(Sig60, (Tag == QueryTag & Genotype == "T"))
    Tag60Susceptible <- subset(Sig60, (Tag == QueryTag & Genotype == "S"))
    
    # Em 90 dias
    Tag90Tolerant <- subset(Sig90, (Tag == QueryTag & Genotype == "T"))
    Tag90Susceptible <- subset(Sig90, (Tag == QueryTag & Genotype == "S"))
    
    # Quais clusters?
    clusterA <- input$CLA
    clusterB <- input$CLB
    
    # Seleção do kMeans para o cluster
    if (InputDate == 30) {
      kmeans <- kMeans30
      SigDay <- Sig30
      PickedTagTolerant <- Tag30Tolerant
      PickedTagSusceptible <- Tag30Susceptible
    } else if (InputDate == 60) {
      kmeans <- kMeans60
      SigDay <- Sig60
      PickedTagTolerant <- Tag60Tolerant
      PickedTagSusceptible <- Tag60Susceptible
    } else if (InputDate == 90) {
      kmeans <- kMeans90
      SigDay <- Sig90
      PickedTagTolerant <- Tag90Tolerant
      PickedTagSusceptible <- Tag90Susceptible
    }
    
    # Montagem da matriz de entrada
    kmeans$centers <- cbind(kmeans$centers, 30, 60, 90)
    colnames(kmeans$centers)[4:6] <- c("I30", "I60", "I90")
    mainDF <- as.data.frame(kmeans$centers)
    DFSubset <- subset(mainDF, rownames(mainDF) == clusterA | rownames(mainDF) == clusterB)
    
    # Separação dos clusters
    clusterA.values <- mainDF[clusterA,1:6]
    clusterB.values <- mainDF[clusterB,1:6]
    
    # Armazenamento das cores (para escala)
    ScaleColors <- c("Cluster A" = "blue",
                     "Cluster B" = "tomato",
                     "Tolerante" = "limegreen",
                     "Suscetivel" = "orangered4")
    
    # Forja do plot básico
    ggplot(mainDF) +
      # Pontos dos clusters a cada período, por cores
      ## Cluster A
      geom_point(aes(x=I30[clusterA], y=R30[clusterA]), colour="blue") +
      geom_point(aes(x=I60[clusterA], y=R60[clusterA]), colour="blue") +
      geom_point(aes(x=I90[clusterA], y=R90[clusterA]), colour="blue") +
      
      ## Cluster B
      geom_point(aes(x=I30[clusterB], y=R30[clusterB]), colour="tomato") +
      geom_point(aes(x=I60[clusterB], y=R60[clusterB]), colour="tomato") +
      geom_point(aes(x=I90[clusterB], y=R90[clusterB]), colour="tomato") +
      
      # Rótulos dos gráficos e estéticas adicionais
      xlab("Tempo (dias)") +
      ylab("Comportamento em relação à não-estressada") +
      geom_hline(mapping=aes(yintercept=0), linetype="dotted") +
      
      # Linha de tendência dos clusters selecionados
      ### Cluster A
      ## 30 pra 60
      geom_line(data=DFSubset,
                mapping=aes(x=c(30, 60),
                            y=c(R30[1], R60[1]),
                            colour="Cluster A"),
                linetype="twodash") +
      
      ## 60 pra 90, cluster A
      geom_line(data=DFSubset,
                mapping=aes(x=c(60, 90),
                            y=c(R60[1], R90[1])),
                colour="blue",
                linetype="twodash") +
      
      ### Cluster B
      ## 30 pra 60
      geom_line(data=DFSubset,
                mapping=aes(x=c(30, 60),
                            y=c(R30[2], R60[2]),
                            colour="Cluster B"),
                linetype="twodash") +
      
      ## 60 pra 90, cluster A
      geom_line(data=DFSubset,
                mapping=aes(x=c(60, 90),
                            y=c(R60[2], R90[2])),
                colour="tomato",
                linetype="twodash") +
      
      ## Verificando valores da tag para suscetível e tolerante
      ## Adicionando os pontos para tag escolhida, tolerante e suscetível
      
      
      ## TOLERANTE ##
      # Tag escolhida em 30 (Tolerant)
      geom_point(data=PickedTagTolerant, mapping=aes(x=30,y=as.double(PickedTagTolerant$R30)),
                 colour=colorPickerA) + 
      
      # Tag escolhida em 60
      geom_point(data=PickedTagTolerant, mapping=aes(x=60,y=as.double(PickedTagTolerant$R60)),
                 colour=colorPickerA) + 
      
      # Tag escolhida em 60
      geom_point(data=PickedTagTolerant, mapping=aes(x=90,y=as.double(PickedTagTolerant$R90)),
                 colour=colorPickerA) + 
      
      # Linha de 30 para 60, tolerante
      geom_line(data=PickedTagTolerant[1:2,],
                mapping=aes(x=c(30, 60),
                            y=c(as.double(PickedTagTolerant$R30), as.double(PickedTagTolerant$R60)),
                            colour="Tolerante"), size=1.5) +
      
      # Linha de 60 para 90, tolerante
      geom_line(data=PickedTagTolerant[2:3,],
                mapping=aes(x=c(60, 90),
                            y=c(as.double(PickedTagTolerant$R60), as.double(PickedTagTolerant$R90)),
                            colour="Tolerante"), size=1.5) +
      
      ## SUSCETÍVEL ##
      # Tag escolhida em 30, Suscetível
      geom_point(data=PickedTagSusceptible, mapping=aes(x=30,y=as.double(PickedTagSusceptible$R30)),
                 colour=colorPickerB) + 
      
      # Tag escolhida em 60
      geom_point(data=PickedTagSusceptible, mapping=aes(x=60,y=as.double(PickedTagSusceptible$R60)),
                 colour=colorPickerB) + 
      
      # Tag escolhida em 60
      geom_point(data=PickedTagSusceptible, mapping=aes(x=90,y=as.double(PickedTagSusceptible$R90)),
                 colour=colorPickerB) + 
      
      # Linha de 30 para 60, Suscetível
      geom_line(data=PickedTagSusceptible[1:2,],
                mapping=aes(x=c(30, 60),
                            y=c(as.double(PickedTagSusceptible$R30), as.double(PickedTagSusceptible$R60)),
                            colour="Suscetivel"), size=1.5) +
      
      # Linha de 60 para 90, tolerante
      geom_line(data=PickedTagSusceptible[2:3,],
                mapping=aes(x=c(60, 90),
                            y=c(as.double(PickedTagSusceptible$R60), as.double(PickedTagSusceptible$R90)),
                            colour="Suscetivel"), size=1.5) +
      
      scale_colour_manual(name="Legenda",
                          values=ScaleColors)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

