## Script completo: Das contagens até o cluster
library("BiocParallel")
library("tximport")
library("DESeq2")
library("clValid")
library("ggplot2")

# (Rafael) Ler:
# https://support.bioconductor.org/p/63201/


# Seções do código
# 1. Entradas do usuário
# 2. Processamento da expressão diferencial (DESeq)
# 3. Realização dos agrupamentos e criação da matriz por genótipos


## [.. 1. Entradas ..] ##

# Objetos:
## dir: pasta onde estão localizados os resultados de quantificação do Salmon
## samples: arquivo de texto que contem o nome das pastas (uma pasta por amostra) e os
# desenhos experimentais
## design -> o que se quer explorar
# efeitos do genótipo apenas: design <- c("genotype")
# interação do genótipo com severidade: c("genotype * severity")
# interação total: design <- c("genotype * is.stressed * time")
## nclust -> número de agrupamentos kMeans
## CutValue: valor de corte em relação a não estressado; já estando diferencialmente expresso;

print("Initializing. Reading sample files and metadata (txt)")
## Entradas:
dir <- "/home/rcsilva/projects/dtr/DESeq2-good"
samples <- read.table("/home/rcsilva/projects/dtr/DESeq2-good/samples.txt", header=T)
nclust = 09
CutValue = 02

# Como consultar interações:
# http://seqanswers.com/forums/showthread.php?t=49527 - M. Love

# Alocar processadores para processos em paralelo, Thor: 30
register(MulticoreParam(20))

## Lendo cada um dos arquivos de quantificação, na pasta dir
files <- file.path(dir, "salmon", samples$folder, "quant.sf")

## Adiciona nomes às amostras
names(files) <- (samples$folder)

## Objeto do tximport para importar no DESeq2
txi <- tximport(files, type = "salmon", txOut=TRUE)

## Adicionando o par binário "estressado/não-estressado"
#samples[,dim(samples)[2]+1] <- NA
samples[,7] <- NA
counter = 1;

for (stress.level in samples$severity) {
  if (stress.level == "Moderate" | stress.level == "Severe") {
    samples[counter,7] <- c("Yes");
  } else {
    samples[counter,7] <- c("No");
  }
  counter = counter+1
}

# Alterando o nome da coluna 07 para is.stressed
names(samples)[7] <- c("is.stressed")
# Tornar um fator:
samples$is.stressed <- as.factor(samples$is.stressed)

# Vinculando genótipo e tempo
samples$genotype_time <- factor(paste0(samples$genotype, "_", samples$time))

# Vinculando genótipo, tempo, e se está estressado
samples$genotype_time_stress <- factor(paste0(samples$genotype_time, "_", samples$is.stressed))



## [.. 2. Processamento DESeq2 ..] ##

## Criando objeto (dds) para o DESeq2 para TODOS os fatores
dds <- DESeqDataSetFromTximport(txi, colData = samples,
                                design =~ genotype_time_stress)

## Colapsar as réplicas técnicas
dds <- collapseReplicates(dds, dds$sample)
                                  
## Rodando o DESeq2
print("Running DESeq2")
print("Choosen contrast:")
print(dds@design)

deseq.dds <- DESeq(dds, parallel=T)

# Inicializa a variável de resultados
results.contrast <- list()

# Armazena os DEs (tolerante vs. suscetível, 30d)
results.contrast$TvS30D <- results(deseq.dds,
                                   # Escolhe um dos contrastes
                                   # com base no resultsNames
                                   contrast=c("genotype_time_stress",
                                              "Tolerant_30d_Yes",
                                              "Susceptible_30d_Yes"),
                                   parallel = T)

# Armazena os DEs (tolerante vs. suscetível, 60d)
results.contrast$TvS60D <- results(deseq.dds,
                                   # Escolhe um dos contrastes
                                   # com base no resultsNames
                                   contrast=c("genotype_time_stress",
                                              "Tolerant_60d_Yes",
                                              "Susceptible_60d_Yes"),
                                   parallel = T)

# Armazena os DEs (tolerante vs. suscetível, 90d)
results.contrast$TvS90D <- results(deseq.dds,
                            # Escolhe um dos contrastes
                            # com base no resultsNames
                            contrast=c("genotype_time_stress",
                                       "Tolerant_90d_Yes",
                                       "Susceptible_90d_Yes"),
                            parallel = T)

# Sumário: 
# summary(results.contrast$TvS90D)
# summary(results.contrast$TvS60D)
# summary(results.contrast$TvS30D)


## Cortar só os significativos (p <0.05) : 
# Inicializa variável
results.sig <- list()
# Para 30D
#results.sig$TvS30D <- subset(results.contrast$TvS30D, (padj < 0.05))

results.sig$TvS30D <- subset(results.contrast$TvS30D,
                             (results.contrast$TvS30D$log2FoldChange > 5) | results.contrast$TvS30D$log2FoldChange < -5)
# Para 60D
results.sig$TvS60D <- subset(results.contrast$TvS60D,
                             (results.contrast$TvS60D$log2FoldChange > 5) | results.contrast$TvS60D$log2FoldChange < -5)
# Para 90D
results.sig$TvS90D <- subset(results.contrast$TvS90D,
                             (results.contrast$TvS90D$log2FoldChange > 5) | results.contrast$TvS90D$log2FoldChange < -5)

## [.. Final do trecho do DESeq2 ..] ##
## [.. Início da obtenção dos clusters ..] ##

print("DESeq2 OK (p < 0.05)")
print("Initializing Clusters/kMeans")

contrast = 1;
ListElement = 1;
# 
clusters <- list();

for (contrast in 1:length(results.sig)) {
  # Teste do contraste:
  print(contrast);
  # Seleção da matriz DE
  print("Capturing significant tags from DE Matrix")
  DESeqResults <- results.sig[[contrast]];
  
  ddfx.gt <- matrix(nrow=length(DESeqResults$baseMean), ncol=3)
  ddfx.gs <- matrix(nrow=length(DESeqResults$baseMean), ncol=3)
  counter = 1;
  
  # Progress messages
  print("Calculating L2FC for each time point")
  IsoformNo <- 1
  
  for (isoform in DESeqResults@rownames){
    
    print("Isoform no:")
    print(IsoformNo)
    print("From a total of:")
    print(length(DESeqResults@rownames))
    
    # Captura dos valores - médias para cada severidade, por tag
    # 30 dias tolerante
    tag30NST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "30d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    tag30SST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "30d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    # 30 dias suscetível
    
    tag30NSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "30d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Susceptible"])
    
    tag30SSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "30d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Susceptible"])
    
    # 60 dias tolerante
    
    tag60NST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "60d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    tag60SST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "60d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    # 60 dias suscetível
    
    tag60NSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "60d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Susceptible"])
    
    tag60SSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "60d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Susceptible"])
    
    # 90 dias tolerante
    
    tag90NST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "90d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    tag90SST <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "90d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Tolerant"])
    
    # 90 dias suscetível
    tag90NSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "90d" &
        deseq.dds@colData$severity == "Not_stressed" &
        deseq.dds@colData$genotype == "Susceptible"])
    
    tag90SSS <- mean(deseq.dds@assays$data$counts[isoform,][
      deseq.dds@colData$time == "90d" &
        deseq.dds@colData$severity == "Severe" &
        deseq.dds@colData$genotype == "Susceptible"])  
    
    ### Fim da obtenção dos valores
    
    ## Cálculo das razões
    # Tolerante a 30 dias
    
    # [.. Fórmula do L2FC utilizada ..]
    # Isoforma A | Isoforma B
    # L2FC(AB): ((A-B+1)/(A+1))-1
    
    # Adicionada pseudocontagem (+1) em ambos os lados para remover
    # a geração de infinitos
    
    T30D <- log2( (tag30SST+1) / (tag30NST+1))
    T60D <- log2( (tag60SST+1) / (tag60NST+1))
    T90D <- log2( (tag90SST+1) / (tag90NST+1))
    S30D <- log2( (tag30SSS+1) / (tag30NSS+1))
    S60D <- log2( (tag60SSS+1) / (tag60NSS+1))
    S90D <- log2( (tag90SSS+1) / (tag90NSS+1))
    
    # Criação da matriz:
    # Adiciona uma linha com valores para a tag na matriz
    ddfx.gt[counter,] <- rbind(T30D, T60D, T90D)
    ddfx.gs[counter,] <- rbind(S30D, S60D, S90D)
    
    # Soma ao contador                        
    counter = counter + 1
    IsoformNo <- IsoformNo + 1
  }
  
  # Adição à lista
  print("Binding L2FC matrix to list structure")
  
  # Alocar espaços 3, 4 da lista
  if (contrast == 2) {
    ListElement <- 3
  }
  
  # Alocar espaços 5, 6 da lista
  if (contrast == 3) {
    ListElement <- 5
  }
  
  clusters[[ListElement]] <- ddfx.gt
  clusters[[ListElement+1]] <- ddfx.gs
  
  # Adicionando nomes às linhas
  print("Adding names to tags...")
  rownames(clusters[[ListElement]]) <- paste0(DESeqResults@rownames, ":", "T")
  rownames(clusters[[ListElement+1]]) <- paste0(DESeqResults@rownames, ":", "S")
  colnames(clusters[[ListElement]]) <- c("R30", "R60", "R90")
  colnames(clusters[[ListElement+1]]) <- c("R30", "R60", "R90")
}

print("L2FC capture ok")
## [.. Final da captura da matriz DE ..] ##

## [.. Organização dos clusters (kMeans) ..] ##
print("Reunindo os clusters em uma matriz única por data...")
## Agrupando os tolerantes e suscetíveis em 30, 60, 90

# Inicializando a lista
AmbosClusters <- list()

## Agrupando DEs por data:
# 1 e 2 para 30 dias
# 3 e 4 para 60 dias
# 5 e 6 para 90 dias

AmbosClusters$at30 <- rbind(clusters[[1]], clusters[[2]])
AmbosClusters$at60 <- rbind(clusters[[3]], clusters[[4]])
AmbosClusters$at90 <- rbind(clusters[[5]], clusters[[6]])

## Teste de filtragem: significativos
# Inicialização
rm(AmbosClustersSig)
AmbosClustersSig <- list();

# Significativos a 30d
AmbosClustersSig$at30 <- AmbosClusters$at30#[AmbosClusters$at30[,1] >= CutValue |
                                            #  AmbosClusters$at30[,1] <= -CutValue, ]

# Significativos a 60d
AmbosClustersSig$at60 <- AmbosClusters$at60#[AmbosClusters$at60[,2] >= CutValue |
                                            #  AmbosClusters$at60[,2] <= -CutValue, ]
# Significativo a 90d 
AmbosClustersSig$at90 <- AmbosClusters$at90#[AmbosClusters$at90[,3] >= CutValue |
                                           #AmbosClusters$at90[,3] <= -CutValue, ]

print("Executando k-means...")
# Montagem de kMeans
kMeans30 <- kmeans(AmbosClustersSig$at30, nclust)
kMeans60 <- kmeans(AmbosClustersSig$at60, nclust)
kMeans90 <- kmeans(AmbosClustersSig$at90, nclust)
print("Atribuindo o número do cluster...")

# Atribuição do cluster na matriz original
AmbosClustersSig$at30 <- cbind(AmbosClustersSig$at30, as.matrix(kMeans30$cluster))
AmbosClustersSig$at60 <- cbind(AmbosClustersSig$at60, as.matrix(kMeans60$cluster))
AmbosClustersSig$at90 <- cbind(AmbosClustersSig$at90, as.matrix(kMeans90$cluster))

# Adicionando tag e genótipo à matriz
rm(genovector)
rm(tags)
rm(genotypes)
# Colando na matriz
for (timet in 1:length(names(AmbosClustersSig))) {
  genovector <- unlist(strsplit(names(AmbosClustersSig[[timet]][,1]), ":")[1:length(AmbosClustersSig[[timet]])])
  tags <- genovector[c(T, F)]
  genotypes <- genovector[c(F, T)]
  AmbosClustersSig[[timet]] <- cbind(AmbosClustersSig[[timet]], as.matrix(tags))
  AmbosClustersSig[[timet]] <- cbind(AmbosClustersSig[[timet]], as.matrix(genotypes))
}

# Adicionando nomes às colunas
colnames(AmbosClustersSig[[1]])[4:6] <- c("Cluster", "Tag", "Genotype")
colnames(AmbosClustersSig[[2]])[4:6] <- c("Cluster", "Tag", "Genotype")
colnames(AmbosClustersSig[[3]])[4:6] <- c("Cluster", "Tag", "Genotype")

  ### [.. Geração de Saídas ..] ###

  # A. Matrizes LFC em relação a não estressado

  # Matriz significante a 30 dias (GT vs GS)
  output.significant.at30 <- as.data.frame(AmbosClustersSig$at30)

  # Matriz significante a 60 dias (GT vs GS)
  output.significant.at60 <- as.data.frame(AmbosClustersSig$at60)

  # Matriz significante a 90 dias (GT vs GS)
  output.significant.at90 <- as.data.frame(AmbosClustersSig$at90)
  
  ## Lista com as tags que não pertencem ao mesmo cluster ##
  # para 30 dias
  output.significant.at30.z <- output.significant.at30[, c('Cluster','Genotype','Tag')]
  
  output.significant.at30.rsp <- reshape(output.significant.at30.z, v.names='Cluster',
                                         idvar='Tag',
                                         timevar='Genotype',
                                         direction="wide")
  
  rownames(output.significant.at30.rsp) <- output.significant.at30.rsp$Tag
  diff.at30 <- subset(output.significant.at30.rsp, output.significant.at30.rsp[,2] != output.significant.at30.rsp[,3] )
  
  # para 60 dias
  output.significant.at60.z <- output.significant.at60[, c('Cluster','Genotype','Tag')]
  
  output.significant.at60.rsp <- reshape(output.significant.at60.z, v.names='Cluster',
                                         idvar='Tag',
                                         timevar='Genotype',
                                         direction="wide")
  
  rownames(output.significant.at60.rsp) <- output.significant.at60.rsp$Tag
  diff.at60 <- subset(output.significant.at60.rsp, output.significant.at60.rsp[,2] != output.significant.at60.rsp[,3] )
  

  # para 90 dias
  output.significant.at90.z <- output.significant.at90[, c('Cluster','Genotype','Tag')]
  
  output.significant.at90.rsp <- reshape(output.significant.at90.z, v.names='Cluster',
                                         idvar='Tag',
                                         timevar='Genotype',
                                         direction="wide")
  
  rownames(output.significant.at90.rsp) <- output.significant.at90.rsp$Tag
  diff.at90 <- subset(output.significant.at90.rsp, output.significant.at90.rsp[,2] != output.significant.at90.rsp[,3] )
  
  ## [.. Lista final de saídas: ..]
  
  # A 30 dias:
  output.significant.at30       # Matriz de clusters para 30 dias
  output.significant.at30.rsp   # Consulta aos clusters por genótipo
  diff.at30                     # Tags que não pertencem ao mesmo cluster
  
  # A 60 dias:
  output.significant.at60       # Matriz de clusters para 30 dias
  output.significant.at60.rsp   # Consulta aos clusters por genótipo
  diff.at60                     # Tags que não pertencem ao mesmo cluster
  
  # A 90 dias:
  output.significant.at90       # Matriz de clusters para 30 dias
  output.significant.at90.rsp   # Consulta aos clusters por genótipo
  diff.at90                     # Tags que não pertencem ao mesmo cluster
