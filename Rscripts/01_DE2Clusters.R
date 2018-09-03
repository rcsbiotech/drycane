## Script completo: Das contagens até o cluster
library("BiocParallel")
library("tximport")
library("DESeq2")
library("clValid")
library("ggplot2")
library("xlsx")

# Code sections
# 1. User inputs
# 2. Differential expression process with DESeq2
# 3. Clustering by genotypes and stress

## [.. 1. Inputs ..] ##

# Objects:
## dir: folder with salmon quantification files
## samples: text file with experimental design and folder names
## (WIP) design -> binary experimental design matrix
## parameter.nclust -> number of desired kMeans clusters
## parameter.cutvalue -> cut-value over not stressed
## parameter.fc.threshold: value to split clustering by signal

print("Initializing. Reading sample files and metadata (txt)")

## Entradas:
dir <- "/home/rcsilva/projects/dtr/DESeq2-good"

# Amostras
samples <- read.table("/home/rcsilva/projects/dtr/DESeq2-good/samples.txt", header=T)

# Corte de FC no DESeq2
parameter.lfc.sig <- 04

# Número de Clusters kMeans
parameter.nclust <- 09

# FC em relação à não estressada para excluir
parameter.cutvalue <- 00

# Valor para separação dos grupos de sinais;
parameter.fc.threshold <- 0.5

# Número de processadores, não usar mais que 20 no Thor
parameter.threads <- 35

### FINALIZAÇÃO DOS INPUTS







### INÍCIO DO SCRIPT

# Como consultar interações:
# http://seqanswers.com/forums/showthread.php?t=49527 - M. Love

# Alocar processadores para processos em paralelo, Thor: 30
register(MulticoreParam(parameter.threads))

## Lendo cada um dos arquivos de quantificação, na pasta dir
files <- file.path(dir, "salmon", samples$folder, "quant.sf")

## Adiciona nomes às amostras
names(files) <- (samples$folder)

## Objeto do tximport para importar no DESeq2
txi <- tximport(files, type = "salmon", txOut=TRUE)

## Adicionando o par binário "estressado/não-estressado"
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
print("Creating DESeq2 object")
dds <- DESeqDataSetFromTximport(txi, colData = samples,
                                design =~ genotype_time_stress)

## Colapsar as réplicas técnicas
print("Collapsing technical replicates")
dds <- collapseReplicates(dds, dds$sample)
                                  
## Rodando o DESeq2
print("Running DESeq2")
print("Choosen contrast:")
print(dds@design)

# Execução do DESeq
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


## Cortar só os significativos (p <0.05) : 
# Inicializa variável
results.sig <- list()
# Para 30D
#results.sig$TvS30D <- subset(results.contrast$TvS30D, (padj < 0.05))

results.sig$TvS30D <- subset(results.contrast$TvS30D,
                             (results.contrast$TvS30D$log2FoldChange > parameter.lfc.sig) |
                               results.contrast$TvS30D$log2FoldChange < -parameter.lfc.sig)
# Para 60D
results.sig$TvS60D <- subset(results.contrast$TvS60D,
                             (results.contrast$TvS60D$log2FoldChange > parameter.lfc.sig) |
                               results.contrast$TvS60D$log2FoldChange < -parameter.lfc.sig)
# Para 90D
results.sig$TvS90D <- subset(results.contrast$TvS90D,
                             (results.contrast$TvS90D$log2FoldChange > 4 ) |
                               results.contrast$TvS90D$log2FoldChange < -4)

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
    # a geração de infinitos positivos e negativos
    
    T30D <- log2( (tag30SST+1) / (tag30NST+1) )
    T60D <- log2( (tag60SST+1) / (tag60NST+1) )
    T90D <- log2( (tag90SST+1) / (tag90NST+1) )
    S30D <- log2( (tag30SSS+1) / (tag30NSS+1) )
    S60D <- log2( (tag60SSS+1) / (tag60NSS+1) )
    S90D <- log2( (tag90SSS+1) / (tag90NSS+1) )
    
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
AmbosClustersSig$at30 <- AmbosClusters$at30[AmbosClusters$at30[,1] >= parameter.cutvalue |
                                              AmbosClusters$at30[,1] <= -parameter.cutvalue, ]

# Significativos a 60d
AmbosClustersSig$at60 <- AmbosClusters$at60[AmbosClusters$at60[,2] >= parameter.cutvalue |
                                              AmbosClusters$at60[,2] <= -parameter.cutvalue, ]
# Significativo a 90d 
AmbosClustersSig$at90 <- AmbosClusters$at90[AmbosClusters$at90[,3] >= parameter.cutvalue |
                                           AmbosClusters$at90[,3] <= -parameter.cutvalue, ]

print("Executando k-means...")
# Montagem de kMeans
kMeans30 <- kmeans(AmbosClustersSig$at30, parameter.nclust, nstart=100)
kMeans60 <- kmeans(AmbosClustersSig$at60, parameter.nclust, nstart=100)
kMeans90 <- kmeans(AmbosClustersSig$at90, parameter.nclust, nstart=100)
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
  
  # "Reforja" a matriz, separando as tags por genótipo por cluster
  output.significant.at30.rsp <- reshape(output.significant.at30.z, v.names='Cluster',
                                         idvar='Tag',
                                         timevar='Genotype',
                                         direction="wide")
  
  # Adiciona nomes às linhas
  rownames(output.significant.at30.rsp) <- output.significant.at30.rsp$Tag
  
  # Separa somente aqueles que são de clusters diferentes
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
  
  # Adiciona nomes às linhas
  rownames(output.significant.at90.rsp) <- output.significant.at90.rsp$Tag
  
  # Isola somente aqueles em clusters diferentes
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
  
  
  
  
  ## [.. Gerando perfis por sinal biológico ..] ##
  
  # Entrada das matrizes
  Sig30 <- as.data.frame(AmbosClustersSig$at30, stringsAsFactors=F)
  Sig60 <- as.data.frame(AmbosClustersSig$at60, stringsAsFactors=F)
  Sig90 <- as.data.frame(AmbosClustersSig$at90, stringsAsFactors=F)
  
  # Montagem dos perfis
  # P - plus (diferença > +0.5 )
  # E - equal (diferença < |0.5| )
  # M - minus (diferença < -0.5 ), assim
  # PP - Sobe em todos os tempos
  # PE - Sobe, e normaliza
  # PM - Sobe, e desce
  # EP - Igual, e sobe
  # EM - Igual, e desce
  # EE - Igual, igual
  # ME - Desce, e equaliza
  # MP - Desce, e sobe
  # MM - Desce, desce
  
  
  
  
  
  ## A 30 dias ###
  ## Comportamento PX (Plus, e algo) ##
  
  print("Criando perfis a 30 dias")
  
  # Curva maior que 0.5 para os dois pontos
  Sig30PP <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig30[,3]) - as.double(Sig30[,2]) > parameter.fc.threshold))))
  # Exportar como excel
  write.table(Sig30PP, file="Sig30PP.txt", sep="\t", quote=F)
  
  
  # Curva maior que meio para o primeiro ponto
  Sig30PE <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) > parameter.fc.threshold)) &
                              (((as.double(Sig30[,3]) - as.double(Sig30[,2]) < parameter.fc.threshold) &
                                  as.double(Sig30[,3]) - as.double(Sig30[,2]) > -parameter.fc.threshold))))
  write.table(Sig30PE, file="Sig30PE.txt", sep="\t", quote=F)
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig30PM <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig30[,3]) - as.double(Sig30[,2]) < -parameter.fc.threshold))))
  write.table(Sig30PM, file="Sig30PM.txt", sep="\t", quote=F)
  

    ## Comportamento EX (equal, e algo)
    # Igual, igual (EE)
  
  Sig30EE <- subset(Sig30, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 30
    (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < parameter.fc.threshold) &
        as.double(Sig30[,2]) - as.double(Sig30[,1]) > -parameter.fc.threshold)))  &
      # Diferenças menores que meio, e maiores que menos meio para 60 e 30
      (((as.double(Sig30[,3]) - as.double(Sig30[,2]) < parameter.fc.threshold) &
          as.double(Sig30[,3]) - as.double(Sig30[,2]) > -parameter.fc.threshold)))
  write.table(Sig30EE, file="Sig30EE.txt", sep="\t", quote=F)
  
  # Igual, e sobe;
  
  Sig30EP <- subset(Sig30, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 30
    (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < parameter.fc.threshold) &
        as.double(Sig30[,2]) - as.double(Sig30[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig30[,3]) - as.double(Sig30[,2]) > parameter.fc.threshold)))
  write.table(Sig30EP, file="Sig30EP.txt", sep="\t", quote=F)
  
  # Igual, e desce;
  
  Sig30EM <- subset(Sig30, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 30
    (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < parameter.fc.threshold) &
        as.double(Sig30[,2]) - as.double(Sig30[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig30[,3]) - as.double(Sig30[,2]) < -parameter.fc.threshold)))
  write.table(Sig30EM, file="Sig30EM.txt", sep="\t", quote=F)
  
  ## Comportamento MX (minus, e algo)
  # Minus, minus (MM)
  # Curva maior que 0.5 para os dois pontos
  Sig30MM <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig30[,3]) - as.double(Sig30[,2]) < -parameter.fc.threshold))))
  write.table(Sig30MM, file="Sig30MM.txt", sep="\t", quote=F)
  
  # Curva maior que meio para o primeiro ponto
  Sig30ME <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < -parameter.fc.threshold)) &
                              (((as.double(Sig30[,3]) - as.double(Sig30[,2]) < parameter.fc.threshold) &
                                  as.double(Sig30[,3]) - as.double(Sig30[,2]) > -parameter.fc.threshold))))
  write.table(Sig30ME, file="Sig30ME.txt", sep="\t", quote=F)
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig30MP <- subset(Sig30, (((as.double(Sig30[,2]) - as.double(Sig30[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig30[,3]) - as.double(Sig30[,2]) > parameter.fc.threshold))))
  write.table(Sig30MP, file="Sig30MP.txt", sep="\t", quote=F)
  
  
  ### Fim: 30 dias ###
  ### A 60 dias ###
  ## Comportamento PX (Plus, e algo) ##
  
  # Curva maior que 0.5 para os dois pontos
  Sig60PP <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig60[,3]) - as.double(Sig60[,2]) > parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro ponto
  Sig60PE <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) > parameter.fc.threshold)) &
                              (((as.double(Sig60[,3]) - as.double(Sig60[,2]) < parameter.fc.threshold) &
                                  as.double(Sig60[,3]) - as.double(Sig60[,2]) > -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig60PM <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig60[,3]) - as.double(Sig60[,2]) < -parameter.fc.threshold))))
  
  
  ## Comportamento EX (equal, e algo)
  # Igual, igual (EE)
  
  Sig60EE <- subset(Sig60, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 60
    (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < parameter.fc.threshold) &
        as.double(Sig60[,2]) - as.double(Sig60[,1]) > -parameter.fc.threshold)))  &
      # Diferenças menores que meio, e maiores que menos meio para 60 e 60
      (((as.double(Sig60[,3]) - as.double(Sig60[,2]) < parameter.fc.threshold) &
          as.double(Sig60[,3]) - as.double(Sig60[,2]) > -parameter.fc.threshold)))
  
  # Igual, e sobe;
  
  Sig60EP <- subset(Sig60, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 60
    (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < parameter.fc.threshold) &
        as.double(Sig60[,2]) - as.double(Sig60[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig60[,3]) - as.double(Sig60[,2]) > parameter.fc.threshold)))
  
  # Igual, e desce;
  
  Sig60EM <- subset(Sig60, (
    # Diferenças menores que meio, e maiores que menos meio para 60 e 60
    (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < parameter.fc.threshold) &
        as.double(Sig60[,2]) - as.double(Sig60[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig60[,3]) - as.double(Sig60[,2]) < -parameter.fc.threshold)))
  
  ## Comportamento MX (minus, e algo)
  # Minus, minus (MM)
  # Curva maior que 0.5 para os dois pontos
  Sig60MM <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig60[,3]) - as.double(Sig60[,2]) < -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro ponto
  Sig60ME <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < -parameter.fc.threshold)) &
                              (((as.double(Sig60[,3]) - as.double(Sig60[,2]) < parameter.fc.threshold) &
                                  as.double(Sig60[,3]) - as.double(Sig60[,2]) > -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig60MP <- subset(Sig60, (((as.double(Sig60[,2]) - as.double(Sig60[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig60[,3]) - as.double(Sig60[,2]) > parameter.fc.threshold))))
  
  
  
  ### A 90 dias ###
  ## Comportamento PX (Plus, e algo) ##
  
  # Curva maior que 0.5 para os dois pontos
  Sig90PP <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig90[,3]) - as.double(Sig90[,2]) > parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro ponto
  Sig90PE <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) > parameter.fc.threshold)) &
                              (((as.double(Sig90[,3]) - as.double(Sig90[,2]) < parameter.fc.threshold) &
                                  as.double(Sig90[,3]) - as.double(Sig90[,2]) > -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig90PM <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) > parameter.fc.threshold)) &
                              ((as.double(Sig90[,3]) - as.double(Sig90[,2]) < -parameter.fc.threshold))))
  
  
  ## Comportamento EX (equal, e algo)
  # Igual, igual (EE)
  
  Sig90EE <- subset(Sig90, (
    # Diferenças menores que meio, e maiores que menos meio para 90 e 90
    (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < parameter.fc.threshold) &
        as.double(Sig90[,2]) - as.double(Sig90[,1]) > -parameter.fc.threshold)))  &
      # Diferenças menores que meio, e maiores que menos meio para 90 e 90
      (((as.double(Sig90[,3]) - as.double(Sig90[,2]) < parameter.fc.threshold) &
          as.double(Sig90[,3]) - as.double(Sig90[,2]) > -parameter.fc.threshold)))
  
  # Igual, e sobe;
  
  Sig90EP <- subset(Sig90, (
    # Diferenças menores que meio, e maiores que menos meio para 90 e 90
    (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < parameter.fc.threshold) &
        as.double(Sig90[,2]) - as.double(Sig90[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig90[,3]) - as.double(Sig90[,2]) > parameter.fc.threshold)))
  
  # Igual, e desce;
  
  Sig90EM <- subset(Sig90, (
    # Diferenças menores que meio, e maiores que menos meio para 90 e 90
    (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < parameter.fc.threshold) &
        as.double(Sig90[,2]) - as.double(Sig90[,1]) > -parameter.fc.threshold)))  &
      ((as.double(Sig90[,3]) - as.double(Sig90[,2]) < -parameter.fc.threshold)))
  
  ## Comportamento MX (minus, e algo)
  # Minus, minus (MM)
  # Curva maior que 0.5 para os dois pontos
  Sig90MM <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig90[,3]) - as.double(Sig90[,2]) < -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro ponto
  Sig90ME <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < -parameter.fc.threshold)) &
                              (((as.double(Sig90[,3]) - as.double(Sig90[,2]) < parameter.fc.threshold) &
                                  as.double(Sig90[,3]) - as.double(Sig90[,2]) > -parameter.fc.threshold))))
  
  # Curva maior que meio para o primeiro, e desce para o segundo
  Sig90MP <- subset(Sig90, (((as.double(Sig90[,2]) - as.double(Sig90[,1]) < -parameter.fc.threshold)) &
                              ((as.double(Sig90[,3]) - as.double(Sig90[,2]) > parameter.fc.threshold))))

  
  # Cola as colunas para entrada no Shiny
  Names1 <- unname(AmbosClustersSig$at30[,5])
  Names2 <- unname(AmbosClustersSig$at60[,5])
  Names3 <- unname(AmbosClustersSig$at90[,5])
  NamesAll <- unique(c(Names1, Names2, Names3))
  NamesAll <- as.data.frame(NamesAll, stringsAsFactors=F)
  
  # Exportando os outros arquivos para o Shiny
  saveRDS(NamesAll, "/home/rcsilva/DrycaneShiny/names.RDS")
  saveRDS(Sig30, "/home/rcsilva/DrycaneShiny/Sig30.RDS")
  saveRDS(Sig60, "/home/rcsilva/DrycaneShiny/Sig60.RDS")
  saveRDS(Sig90, "/home/rcsilva/DrycaneShiny/Sig90.RDS")
  saveRDS(kMeans30, "/home/rcsilva/DrycaneShiny/kMeans30.RDS")
  saveRDS(kMeans60, "/home/rcsilva/DrycaneShiny/kMeans60.RDS")
  saveRDS(kMeans90, "/home/rcsilva/DrycaneShiny/kMeans90.RDS")
  saveRDS(AmbosClustersSig, "/home/rcsilva/DrycaneShiny/AmbosClustersSig.RDS")
  #AmbosClustersSig<- readRDS("./AmbosClustersSig.RDS")

  
