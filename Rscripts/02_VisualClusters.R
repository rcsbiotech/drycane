# [ .. 1. Shiny User Inputs Processing ..] #

# Limpar arquivos anteriores
rm(kmeans)
rm(mainDF)
rm(DFSubset)
rm(clusterA.values)
rm(clusterB.values)
rm(SigDay)
rm(c(Tag30Tolerant, Tag60Tolerant, Tag90Tolerant,
     Tag30Susceptible, Tag60Susceptible, Tag90Susceptible))

# Input01: Qual data DEG (30, 60, 90)
# Diferencialmente expresso na data..
InputDate <- 30

# Quais clusters?
clusterA <- 07
clusterB <- 09

# Qual tag consultar?
QueryTag = "TRINITY_DN32024_c0_g1_i3"

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

# Escolha de cores (T, S)
colorPickerA = "limegreen"
colorPickerB = "orangered4"

# Separação dos clusters
clusterA.values <- mainDF[clusterA,1:6]
clusterB.values <- mainDF[clusterB,1:6]

# Forja do plot básico
ggplot(mainDF) +
  # Pontos dos clusters a cada período, por cores
  ## Cluster A
  geom_point(aes(x=I30[clusterA], y=R30[clusterA]), color="blue") +
  geom_point(aes(x=I60[clusterA], y=R60[clusterA]), color="blue") +
  geom_point(aes(x=I90[clusterA], y=R90[clusterA]), color="blue") +
  
  ## Cluster B
  geom_point(aes(x=I30[clusterB], y=R30[clusterB]), color="tomato") +
  geom_point(aes(x=I60[clusterB], y=R60[clusterB]), color="tomato") +
  geom_point(aes(x=I90[clusterB], y=R90[clusterB]), color="tomato") +
  
  # Rótulos dos gráficos e estéticas adicionais
  xlab("Tempo (dias)") +
  ylab("Comportamento em relação à não-estressada") +
  geom_hline(mapping=aes(yintercept=0), linetype="dotted") +
  
  # Linha de tendência dos clusters selecionados
  ### Cluster A
  ## 30 pra 60
  geom_line(data=DFSubset,
    mapping=aes(x=c(30, 60),
                y=c(R30[1], R60[1])),
    color="blue",
    linetype="twodash") +

  ## 60 pra 90, cluster A
  geom_line(data=DFSubset,
          mapping=aes(x=c(60, 90),
                      y=c(R60[1], R90[1])),
          color="blue",
          linetype="twodash") +

  ### Cluster B
  ## 30 pra 60
  geom_line(data=DFSubset,
          mapping=aes(x=c(30, 60),
                      y=c(R30[2], R60[2])),
          color="tomato",
          linetype="twodash") +
  
  ## 60 pra 90, cluster A
  geom_line(data=DFSubset,
            mapping=aes(x=c(60, 90),
                        y=c(R60[2], R90[2])),
            color="tomato",
            linetype="twodash") +

  ## Verificando valores da tag para suscetível e tolerante
  ## Adicionando os pontos para tag escolhida, tolerante e suscetível
  
  
  ## TOLERANTE ##
  # Tag escolhida em 30 (Tolerant)
  geom_point(data=PickedTagTolerant, mapping=aes(x=30,y=as.double(PickedTagTolerant$R30)),
             color=colorPickerA) + 

  # Tag escolhida em 60
  geom_point(data=PickedTagTolerant, mapping=aes(x=60,y=as.double(PickedTagTolerant$R60)),
           color=colorPickerA) + 

  # Tag escolhida em 60
  geom_point(data=PickedTagTolerant, mapping=aes(x=90,y=as.double(PickedTagTolerant$R90)),
           color=colorPickerA) + 

  # Linha de 30 para 60, tolerante
  geom_line(data=PickedTagTolerant[1:2,],
            mapping=aes(x=c(30, 60),
                        y=c(as.double(PickedTagTolerant$R30), as.double(PickedTagTolerant$R60))),
            color=colorPickerA, size=2) +

  # Linha de 60 para 90, tolerante
  geom_line(data=PickedTagTolerant[2:3,],
          mapping=aes(x=c(60, 90),
                      y=c(as.double(PickedTagTolerant$R60), as.double(PickedTagTolerant$R90))),
          color=colorPickerA, size=2) +
  
  ## SUSCETÍVEL ##
  # Tag escolhida em 30 (Tolerant)
  geom_point(data=PickedTagSusceptible, mapping=aes(x=30,y=as.double(PickedTagSusceptible$R30)),
           color=colorPickerB) + 
  
  # Tag escolhida em 60
  geom_point(data=PickedTagSusceptible, mapping=aes(x=60,y=as.double(PickedTagSusceptible$R60)),
             color=colorPickerB) + 
  
  # Tag escolhida em 60
  geom_point(data=PickedTagSusceptible, mapping=aes(x=90,y=as.double(PickedTagSusceptible$R90)),
             color=colorPickerB) + 
  
  # Linha de 30 para 60, tolerante
  geom_line(data=PickedTagSusceptible[1:2,],
            mapping=aes(x=c(30, 60),
                        y=c(as.double(PickedTagSusceptible$R30), as.double(PickedTagSusceptible$R60))),
            color=colorPickerB, size=2) +
  
  # Linha de 60 para 90, tolerante
  geom_line(data=PickedTagSusceptible[2:3,],
            mapping=aes(x=c(60, 90),
                        y=c(as.double(PickedTagSusceptible$R60), as.double(PickedTagSusceptible$R90))),
            color=colorPickerB, size=2) +
  
  scale_color_identity(name = 'the fill', guide = 'legend', labels = c('m1'))
