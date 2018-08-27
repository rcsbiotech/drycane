
## Plot per cluster ##

# Sections:
# 1. Shiny inputs
# 2. Plot generation

# [ .. 1. Shiny User Inputs Processing ..] #

# Input01: Qual data DEG (30, 60, 90)
# Diferencialmente expresso na data..
InputDate <- 30

# Quais clusters?
ClusterA <- 01
ClusterB <- 08

# Limites para o eixo Y
Yupper <- 3
Ylower <- -3

# Init modelos lineares
init = T

# Seleção do kMeans para o cluster
if (InputDate == 30) {
  kmeans <- kMeans30
} else if (InputDate == 60) {
  kmeans <- kMeans60
} else if (InputDate == 90) {
  kmeans <- kMeans90
}

kmeans$centers <- cbind(kmeans$centers, 30, 60, 90)
colnames(kmeans$centers)[4:6] <- c("I30", "I60", "I90")

## <.. End of Shiny Input/Processing ..> ##

## [.. 2. Figure generation .. ] ##
ggplot(AmbosClustersSig[[InputDate/30]][,InputDate/30],
       aes(x=c(30),
           y=AmbosClustersSig[[InputDate/30]][,InputDate/30][names(kmeans$cluster == ClusterA |
                                                                        kmeans$cluster == ClusterB)])) +
  
  # Coordenadas cartesianas: ampliação da figura
  coord_cartesian(ylim=c(Ylower, Yupper)) +
  
  # Adição de rótulos
  ylab("Log2FC") +
  xlab("Course of time") +
  
  # Remover o eixo ys
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  
  # Escala do eixo y
  #scale_y_discrete(breaks=c(2), limits=c(-4, 4)) +
  
                ## [.. Adição dos agrupamentos ..] ##
  
  # Adição do primeiro cluster.
  # mapping = de onde sairão os dados;
  # x = c(dia), o valor para o eixo x;
  # y = os dados a serem retirados do clusterA

  geom_point(mapping=
                aes(x=c(30),
                    y=(AmbosClustersSig[[InputDate/30]][,1][names(kmeans$cluster[kmeans$cluster == ClusterA])])),
              size=0.5, #width = 3,
             color="steelblue") +
  
  geom_point(mapping=
                aes(x=c(60),
                    y=(AmbosClustersSig[[InputDate/30]][,2][names(kmeans$cluster[kmeans$cluster == ClusterA])])),
              size=0.5, #width = 3,
             color="steelblue") +
  
  geom_point(mapping=
                aes(x=c(90),
                    y=(AmbosClustersSig[[InputDate/30]][,3][names(kmeans$cluster[kmeans$cluster == ClusterA])])),
              size=0.5, #width = 3,
             color="steelblue") +
  
  # Adição do segundo cluster (mesma lógica)
  
  geom_point(mapping= 
                aes(x=c(30),
                    y=(AmbosClustersSig[[InputDate/30]][,1][names(kmeans$cluster[kmeans$cluster == ClusterB])])),
              size=0.5,
             #width = 3,
             color="tomato3") +
  
  geom_point(mapping=
                aes(x=c(60),
                    y=(AmbosClustersSig[[InputDate/30]][,2][names(kmeans$cluster[kmeans$cluster == ClusterB])])),
              size=0.5,
             #width = 3,
             color="tomato3") +
  
  geom_point(mapping=
                aes(x=c(90),
                    y=(AmbosClustersSig[[InputDate/30]][,3][names(kmeans$cluster[kmeans$cluster == ClusterB])])),
              size=0.5, #width = 3,
             color="tomato3") +
  
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  
                  ## [.. Adição de centróides ..] ##
  
  # Adição do centróide do cluster A
  
  geom_point(mapping=
               aes(x=c(30, 60, 90),
                   y=c(kmeans$centers[ClusterA, 1],
                       kmeans$centers[ClusterA, 2],
                       kmeans$centers[ClusterA, 3])), size=10, color="steelblue") +
  
  # Adição do centróide do cluster B
  geom_point(mapping=
               aes(x=c(30, 60, 90),
                   y=c(kmeans$centers[ClusterB, 1],
                       kmeans$centers[ClusterB, 2],
                       kmeans$centers[ClusterB, 3])), size=10, color="tomato3") +
  
  # Adição da linha
    geom_line(mapping=
               aes(x=c(30, 60, 90),
                   y=c(kmeans$centers[ClusterA, 1],
                       kmeans$centers[ClusterA, 2],
                       kmeans$centers[ClusterA, 3])), size=1, color="steelblue") +
    

    geom_line(mapping=
             aes(x=c(30, 60, 90),
                 y=c(kmeans$centers[ClusterB, 1],
                     kmeans$centers[ClusterB, 2],
                     kmeans$centers[ClusterB, 3])), size=1, color="tomato3")

  # Plotar junto do gráfico, as tags por cluster

  geom_point(mapping=
               aes(x=c(30),
                   y=(AmbosClustersSig[[InputDate/30]][,1][as.data.frame(
                     AmbosClustersSig$at30,stringsAsFactors=F)$Tag == as.data.frame(
                       diff.at30$Tag,stringsAsFactors=F)),
             size=0.5, #width = 3,
             color="steelblue"))


  
