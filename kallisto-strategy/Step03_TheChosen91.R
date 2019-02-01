# Cria os 88 perfis do Cana-Seca, com base nas anotações:

### Bibliotecas
# Adiciona zero a esquerda
library(stringr)

### Funções próprias
# Negação rápida
'%!in%' <- function(x,y)!('%in%'(x,y))

# Comparação de tamanho
'%sizecheck%' <- function(x,y)(length(x) == length(y))

# x.test = c("2", "banana")
# y.test = "potato"
# z.test = c("monkey", "wrench")

# x.test %sizecheck% y.test
# x.test %sizecheck% z.test
### Fim das declarações de funções


## Camada 01: Isoformas DEs estressadas
##    - Estressada contra não-estressada em relação a 30 dias
##    - ~~~~ a 60 dias
##    - ~~~~ a 90 dias

## Camada 02: Dinâmica temporal das isoformas
##    - Pra cada isoforma acima...
  #   > Ela é diferencialmente expressa entre 60 e 30 dias...
    #   > Para a planta estressada?
    #   > Para a planta não-estressada?

## Requer os scripts já rodados
  # Step01DEGsToClusters
  # Step02VersusBlank


###### Etapa 01. DEs estressadas em relação às não estressadas ######

# Os objetos já existem, se os scripts anteriores foram rodados:

# SS3vSN3
# SS6vSN6
# SS9vSN9
# TS3vTN3
# TS6vTN6
# TS9vTN9


###### Etapa 02. DEs cada uma em relação a seu tempo ######
# Neste ponto, é necessário coletar cada um deles nos resultados do DESeq2:

## Na ordem, para as tolerantes
# - TS6 / TS3
# - TS9 / TS3
# - TN6 / TN3
# - TN9 / TN6

results.part2.contrast$TS6vTS3 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Tolerant_Severe_60d",
                                                     "Tolerant_Severe_30d"),
                                          parallel = T)

results.part2.contrast$TS9vTS6 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Tolerant_Severe_90d",
                                                     "Tolerant_Severe_60d"),
                                          parallel = T)

results.part2.contrast$TN6vTN3 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Tolerant_Not_stressed_60d",
                                                     "Tolerant_Not_stressed_30d"),
                                          parallel = T)

results.part2.contrast$TN9vTN6 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Tolerant_Not_stressed_90d",
                                                     "Tolerant_Not_stressed_60d"),
                                          parallel = T)


## E então para as suscetíveis
# - SS6 / SS3
# - SS9 / SS3
# - SN6 / SN3
# - SN9 / SN6

results.part2.contrast$SS6vSS3 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Susceptible_Severe_60d",
                                                     "Susceptible_Severe_30d"),
                                          parallel = T)

results.part2.contrast$SS9vSS6 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Susceptible_Severe_90d",
                                                     "Susceptible_Severe_60d"),
                                          parallel = T)

results.part2.contrast$SN6vSN3 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Susceptible_Not_stressed_60d",
                                                     "Susceptible_Not_stressed_30d"),
                                          parallel = T)

results.part2.contrast$SN9vSN6 <- results(deseq.dds.part2,
                                          contrast=c("gt_sl_time",
                                                     "Susceptible_Not_stressed_90d",
                                                     "Susceptible_Not_stressed_60d"),
                                          parallel = T)

# Então armazenamos os resultados significativos, para tolerantes e suscetíveis

TS6vTS3 <- results.part2.contrast$TS6vTS3 <- subset(
  x=results.part2.contrast$TS6vTS3,subset=(
    results.part2.contrast$TS6vTS3$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$TS6vTS3$pvalue < parameter.deseq.padj |
      results.part2.contrast$TS6vTS3$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$TS6vTS3$pvalue < parameter.deseq.padj))

TS9vTS6 <- results.part2.contrast$TS9vTS6 <- subset(
  x=results.part2.contrast$TS9vTS6,subset=(
    results.part2.contrast$TS9vTS6$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$TS9vTS6$pvalue < parameter.deseq.padj |
      results.part2.contrast$TS9vTS6$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$TS9vTS6$pvalue < parameter.deseq.padj))

TN9vTN6 <- results.part2.contrast$TN9vTN6 <- subset(
  x=results.part2.contrast$TN9vTN6,subset=(
    results.part2.contrast$TN9vTN6$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$TN9vTN6$pvalue < parameter.deseq.padj |
      results.part2.contrast$TN9vTN6$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$TN9vTN6$pvalue < parameter.deseq.padj))

TN6vTN3 <- results.part2.contrast$TN6vTN3 <- subset(
  x=results.part2.contrast$TN6vTN3,subset=(
    results.part2.contrast$TN6vTN3$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$TN6vTN3$pvalue < parameter.deseq.padj |
      results.part2.contrast$TN6vTN3$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$TN6vTN3$pvalue < parameter.deseq.padj))

## Suscetíveis

SS6vSS3 <- results.part2.contrast$SS6vSS3 <- subset(
  x=results.part2.contrast$SS6vSS3,subset=(
    results.part2.contrast$SS6vSS3$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$SS6vSS3$pvalue < parameter.deseq.padj |
      results.part2.contrast$SS6vSS3$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$SS6vSS3$pvalue < parameter.deseq.padj))

SS9vSS6 <- results.part2.contrast$SS9vSS6 <- subset(
  x=results.part2.contrast$SS9vSS6,subset=(
    results.part2.contrast$SS9vSS6$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$SS9vSS6$pvalue < parameter.deseq.padj |
      results.part2.contrast$SS9vSS6$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$SS9vSS6$pvalue < parameter.deseq.padj))

SN9vSN6 <- results.part2.contrast$SN9vSN6 <- subset(
  x=results.part2.contrast$SN9vSN6,subset=(
    results.part2.contrast$SN9vSN6$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$SN9vSN6$pvalue < parameter.deseq.padj |
      results.part2.contrast$SN9vSN6$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$SN9vSN6$pvalue < parameter.deseq.padj))

SN6vSN3 <- results.part2.contrast$SN6vSN3 <- subset(
  x=results.part2.contrast$SN6vSN3,subset=(
    results.part2.contrast$SN6vSN3$log2FoldChange > parameter.lfc.sig &
      results.part2.contrast$SN6vSN3$pvalue < parameter.deseq.padj |
      results.part2.contrast$SN6vSN3$log2FoldChange < -parameter.lfc.sig &
      results.part2.contrast$SN6vSN3$pvalue < parameter.deseq.padj))

# Ficam, assim, todos os conjuntos prontos

## Camada de informação 01: plantas em relação aos seus controles
# SS3vSN3
# SS6vSN6
# SS9vSN9
# TS3vTN3
# TS6vTN6
# TS9vTN9

## Camada de informação 02: cada planta em cada tratamento em relação a sua etapa anterior no tempo
# TS6vTS3
# TS9vTS6
# TN9vTN6
# TN6vTN3
# SS6vSS3
# SS9vSS6
# SN9vSN6
# SN6vSN3

###### Etapa 03. Associação dos perfis ######

# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # # 
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #
# # !!!!!!!!!!!!!!!!!!!!!!!!AVISO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # #

## Iniciando-se somente da planta tolerante, tentando neutralizar o nome 
## das variáveis, para depois trocá-las pela suscetível se houver interesse.

#### Passo 1. Criando a camada 01. ####

## Isso será feito em duas etapas:
# Primeiro: criar a primeira camada de informação, nos 8 perfis:
# Camada 1, perfil 1: diferencialmente expressos nos 3 tempos
# C1_P1: 1 1 1

## Para cada valor positivo, são necessárias duas colunas com valores de DE:
C1_P1_c1 <- subset(
  x=TS3vTN3[c(2,3)],
  as.vector(TS3vTN3@rownames) %in% as.vector(TS6vTN6@rownames) &
    as.vector(TS3vTN3@rownames) %in% as.vector(TS9vTN9@rownames))

C1_P1_c2 <- subset(
  x=TS6vTN6[c(2,3)],
  as.vector(TS6vTN6@rownames) %in% as.vector(TS9vTN9@rownames) &
    as.vector(TS6vTN6@rownames) %in% as.vector(TS3vTN3@rownames))

C1_P1_c3 <- subset(
  x=TS9vTN9[c(2,3)],
  as.vector(TS9vTN9@rownames) %in% as.vector(TS6vTN6@rownames) &
    as.vector(TS9vTN9@rownames) %in% as.vector(TS3vTN3@rownames))

## Cola os valores
C1_P1 <- cbind(c(C1_P1_c1,
                 C1_P1_c2,
                 C1_P1_c3))

colnames(C1_P1) <- c(
  "LFC TS3/TN3",
  "pval TS3/TN3",
  "LFC TS6/TN6",
  "pval TS6/TN6",
  "LFC TS9/TN9",
  "pval TS9/TN9")

# Camada 1, perfil 2: diferencialmente expresso em 30, 60
# C1_P2: 1 1 0
C1_P2_c1 <- subset(x=TS3vTN3[c(2,3)],
                   as.vector(TS3vTN3@rownames) %in% as.vector(TS6vTN6@rownames) &
                     as.vector(TS3vTN3@rownames) %!in% as.vector(TS9vTN9@rownames))

C1_P2_c2 <- subset(x=TS6vTN6[c(2,3)],
                as.vector(TS6vTN6@rownames) %in% as.vector(TS3vTN3@rownames) &
                  as.vector(TS6vTN6@rownames) %!in% as.vector(TS9vTN9@rownames))

C1_P2 <- cbind(c(
  C1_P2_c1,
  C1_P2_c2
))

colnames(C1_P2) <- c(
  "LFC TS3/TN3",
  "pval TS3/TN3",
  "LFC TS6/TN6",
  "pval TS6/TN6")

# Camada 1, perfil 3: diferencialmente expresso em 30, 90
# C1_P3: 1 0 1
C1_P3_c1 <- subset(
  x=TS3vTN3[c(2,3)],
  as.vector(TS3vTN3@rownames) %in% as.vector(TS9vTN9@rownames) &
    as.vector(TS3vTN3@rownames) %!in% as.vector(TS6vTN6@rownames))

C1_P3_c2 <- subset(
  x=TS9vTN9[c(2,3)],
  as.vector(TS9vTN9@rownames) %in% as.vector(TS3vTN3@rownames) &
    as.vector(TS9vTN9@rownames) %!in% as.vector(TS6vTN6@rownames))

C1_P3 <- cbind(c(
  C1_P3_c1,
  C1_P3_c2
))

colnames(C1_P3) <- c(
  "LFC TS3/TN3",
  "pval TS3/TN3",
  "LFC TS9/TN9",
  "pval TS9/TN9")

# C1_P4: 0 1 1
C1_P4_c1 <- subset(
  x=TS6vTN6[c(2,3)],
  as.vector(TS6vTN6@rownames) %in% as.vector(TS9vTN9@rownames) &
    as.vector(TS6vTN6@rownames) %!in% as.vector(TS3vTN3@rownames))

C1_P4_c2 <- subset(
  x=TS9vTN9[c(2,3)],
  as.vector(TS9vTN9@rownames) %in% as.vector(TS6vTN6@rownames) &
    as.vector(TS9vTN9@rownames) %!in% as.vector(TS3vTN3@rownames))

C1_P4 <- cbind(c(
  C1_P4_c1,
  C1_P4_c2
))

colnames(C1_P4) <- c(
  "LFC TS6/TN6",
  "pval TS6/TN6",
  "LFC TS9/TN9",
  "pval TS9/TN9")

# C1_P5: 1 0 0
C1_P5 <- subset(x=TS3vTN3[c(2,3)],
                as.vector(TS3vTN3@rownames)  %!in% as.vector(TS9vTN9@rownames) &
                  as.vector(TS3vTN3@rownames) %!in% as.vector(TS6vTN6@rownames))

colnames(C1_P5) <- c(
  "LFC TS3/TN3",
  "pval TS3/TN3")

# C1_P7: 0 1 0
C1_P6 <- subset(x=TS6vTN6[c(2,3)],
                as.vector(TS6vTN6@rownames)  %!in% as.vector(TS3vTN3@rownames) &
                  as.vector(TS6vTN6@rownames) %!in% as.vector(TS9vTN9@rownames))

colnames(C1_P6) <- c(
  "LFC TS6/TN6",
  "pval TS6/TN6")

# C1_P6: 0 0 1
C1_P7 <- subset(x=TS9vTN9[c(2,3)],
                as.vector(TS9vTN9@rownames)  %!in% as.vector(TS3vTN3@rownames) &
                  as.vector(TS9vTN9@rownames) %!in% as.vector(TS6vTN6@rownames))

colnames(C1_P7) <- c(
  "LFC TS9/TN9",
  "pval TS9/TN9")

# Vamos armazená-los todos em estrutura de lista, para fácil acesso
## e interação em laços
layer1 <- list()

layer1[1] <- C1_P1
layer1[2] <- C1_P2
layer1[3] <- C1_P3
layer1[4] <- C1_P4
layer1[5] <- C1_P5
layer1[6] <- C1_P6
layer1[7] <- C1_P7

## Para acessar cada um dos elementos (ex, camada 5):
# layers[5][[1]]

#### Passo 2. Criando a camada 02. ####
# Agora, criaremos os 13 perfis temporais, com base nos objetos.

# TS6vTS3
# TS9vTS6
# TN9vTN6
# TN6vTN3
# SS6vSS3
# SS9vSS6
# SN9vSN6
# SN6vSN3

# São 13 camadas previstas:
# 1 1 - 4 valores DE, exclui 0
# 1 0 - 2 valores DE, exclui 2
# 0 1 - 2 valores DE, exclui 2
# 1 E - 3 valores DE, exclui 1
# E 1 - 3 valores DE, exclui 1
# E E - 2 valores DE, exclui 2 and so forth...
# 1 C - 3 valores DE
# C 1 - 3 valores DE
# C C - 2 valores DE
# 0 E - 1 valor DE
# E 0 - 1 valor DE
# 0 C - 1 valor DE
# C 0 - 1 valor DE

## Para construí-las com o valor correto, é necessário que partam da matriz
## De expressão diferencial para cada um de seus tempos, e que também sempre
## sejam removidos os valores os quais não se tem interesse na sobreposição.
# 1 1 - 4 valores DE
## Captura os 4 valores

C2_P1_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN6vTN3@rownames))

C2_P1_c2 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN6vTN3@rownames))

C2_P1_c3 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TN6vTN3@rownames))

C2_P1_c4 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN6vTN3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN6vTN3@rownames) %in% as.vector(TN9vTN6@rownames))

## Une as colunas

C2_P1 <- cbind(
  c(
    C2_P1_c1,
    C2_P1_c2,
    C2_P1_c3,
    C2_P1_c4)
)

## Renomeia os passos
colnames(C2_P1) <- c(
  "L2FC S 60/30",
  "pval S 60/30",
  "L2FC S 90/60",
  "pval S 90/60",
  "L2FC N 60/30",
  "pval N 60/30",
  "L2FC N 90/60",
  "pval N 90/60"
)

# 1 0 - 2 valores DE
C2_P2_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %!in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN9vTN6@rownames))

C2_P2_c2 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %!in% as.vector(TS9vTS6@rownames) &
    as.vector(TN6vTN3@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TN9vTN6@rownames))

C2_P2 <- cbind(
  C2_P2_c1,
  C2_P2_c2
)

colnames(C2_P2) <- c(
    "L2FC S 60/30",
    "pval S 60/30",
    "L2FC N 60/30",
    "pval N 60/30"
)

# 0 1 - 2 valores DE
C2_P3_c1 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %!in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P3_c2 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %!in% as.vector(TS6vTS3@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P3 <- cbind(
  C2_P3_c1,
  C2_P3_c2
)

colnames(C2_P3) <- c(
  "L2FC S 90/60",
  "pval S 90/60",
  "L2FC N 90/60",
  "pval N 90/60"
)

# 1 E - 3 valores DE
C2_P4_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN9vTN6@rownames))

C2_P4_c2 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN9vTN6@rownames))

C2_P4_c3 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN6vTN3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TN9vTN6@rownames))

C2_P4 <- cbind(
  C2_P4_c1,
  C2_P4_c2,
  C2_P4_c3
)

colnames(C2_P4) <- c(
  "L2FC S 60/30",
  "pval S 60/30",
  "L2FC S 90/30",
  "pval S 90/30",
  "L2FC N 60/30",
  "pval N 60/30"
)

# C2_P5: E 1 - 3 valores DE
C2_P5_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P5_c2 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P5_c3 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P5 <- cbind(
  C2_P5_c1,
  C2_P5_c2,
  C2_P5_c3
)

colnames(C2_P5) <- c(
  "L2FC S 60/30",
  "pval S 60/30",
  "L2FC S 90/60",
  "pval S 90/60",
  "L2FC N 90/60",
  "pval N 90/60"
)

# C2_P6: E E - 2 valores DE
C2_P6_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN9vTN6@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P6_c2 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P6 <- cbind(
  C2_P6_c1,
  C2_P6_c2
)

colnames(C2_P6) <- c(
  "L2FC S 60/30",
  "pval S 60/30",
  "L2FC S 90/60",
  "pval S 90/60"
)

# C2_P7: 1 C - 3 valores DE
C2_P7_c1 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TS6vTS3@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P7_c2 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P7_c3 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %in% as.vector(TS6vTS3@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TN6vTN3@rownames))

C2_P7 <- cbind(
  C2_P7_c1,
  C2_P7_c2,
  C2_P7_c3
)

colnames(C2_P7) <- c(
  "L2FC S 60/30",
  "pval S 60/30",
  "L2FC S 90/30",
  "pval S 90/30",
  "L2FC N 90/60",
  "pval N 90/60"
)

# C2_P8: C 1 - 3 valores DE
C2_P8_c1 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TS9vTS6@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TS6vTS3@rownames))

C2_P8_c2 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN6vTN3@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TS6vTS3@rownames))

C2_P8_c3 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TN9vTN6@rownames) %in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TS6vTS3@rownames))

C2_P8 <- cbind(
  C2_P8_c1,
  C2_P8_c2,
  C2_P8_c3
)

colnames(C2_P8) <- c(
  "L2FC S 90/60",
  "pval S 90/60",
  "L2FC N 90/30",
  "pval N 90/30",
  "L2FC N 60/30",
  "pval N 60/30"
)

# C2_P9: C C - 2 valores DE
C2_P9_c1 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %in% as.vector(TN9vTN6@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TS9vTS6@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TS6vTS3@rownames))

C2_P9_c2 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %in% as.vector(TN6vTN3@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TS9vTS6@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TS6vTS3@rownames))

C2_P9 <- cbind(
  C2_P9_c1,
  C2_P9_c2
)

colnames(C2_P9) <- c(
  "L2FC N 60/30",
  "pval N 60/30",
  "L2FC N 90/60",
  "pval N 90/60"
)

# C2_P10: 0 E - 1 valor DE
C2_P10 <- subset(
  x=TS9vTS6[c(2,5)],
  as.vector(TS9vTS6@rownames) %!in% as.vector(TN6vTN3@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TN9vTN6@rownames) &
    as.vector(TS9vTS6@rownames) %!in% as.vector(TS6vTS3@rownames))

colnames(C2_P10) <- c(
  "L2FC N 90/60",
  "pval N 90/60"
)

# C2_P11: E 0 - 1 valor DE
C2_P11 <- subset(
  x=TS6vTS3[c(2,5)],
  as.vector(TS6vTS3@rownames) %!in% as.vector(TN6vTN3@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TN9vTN6@rownames) &
    as.vector(TS6vTS3@rownames) %!in% as.vector(TS9vTS6@rownames))

colnames(C2_P11) <- c(
  "L2FC N 60/30",
  "pval N 60/30"
)

# C2_P12: 0 C - 1 valor DE
C2_P12 <- subset(
  x=TN9vTN6[c(2,5)],
  as.vector(TN9vTN6@rownames) %!in% as.vector(TN6vTN3@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TS6vTS3@rownames) &
    as.vector(TN9vTN6@rownames) %!in% as.vector(TS9vTS6@rownames))

colnames(C2_P12) <- c(
  "L2FC N 90/60",
  "pval N 90/60"
)

# C2_P13: C 0 - 1 valor DE
C2_P13 <- subset(
  x=TN6vTN3[c(2,5)],
  as.vector(TN6vTN3@rownames) %!in% as.vector(TN9vTN6@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TS6vTS3@rownames) &
    as.vector(TN6vTN3@rownames) %!in% as.vector(TS9vTS6@rownames))

colnames(C2_P13) <- c(
  "L2FC N 90/60",
  "pval N 90/60"
)

### E então, colamos juntos também em estrutura de lista ordenada,
## para fácil interação com laços:
layer2 <- list()
layer2[1] <- C2_P1
layer2[2] <- C2_P2
layer2[3] <- C2_P3
layer2[4] <- C2_P4
layer2[5] <- C2_P5
layer2[6] <- C2_P6
layer2[7] <- C2_P7
layer2[8] <- C2_P8
layer2[9] <- C2_P9
layer2[10] <- C2_P10
layer2[11] <- C2_P11
layer2[12] <- C2_P12
layer2[13] <- C2_P13

### ERROR CONTROL: Verificação de bugs e sobreposições
cnt_o <- 1
cnt_i <- 2

for (prof in length(rownames(layer2[cnt_o][[1]]))) {
  for (prof.c in length(rownames(layer2[cnt_i][[1]]))) {
    if (prof %sizecheck% prof.c) {
      print("A profile is matching in size, might be an error.")
      print(c("Criminals: ", cnt_o, cnt_i))
    }
    cnt_i <- cnt_i + 1
  }
  cnt_o <- cnt_o + 1
}

cnt_o <- 1
cnt_i <- 2

#### Passo 3. Mesclando as duas camadas de perfis de maneira combinatória ####
# Agora com as duas camadas prontas em lista, podemos combiná-las.

# As listas são:
## Camada 01
# layer1
## Camada 02
# layer2

# Criaremos um novo objeto em lista (profiles) que irá armazenar
# as combinações únicas de tags entre os dois perfis.

# A lógica do objeto básico, sem estrutura de loop é da seguinte forma:
# Todas as luzes acesas (1 1 1 1 1)

# Primeiro, descobre-se o que eles tem em comum:
test.common <- row.names(
  subset(
    x=layer1[1][[1]],
    row.names(layer1[1][[1]]) %in% row.names(layer2[1][[1]])
  )
)

# Então, usando este nome, cria o perfil para as isoformas em comum
profiles.test <- list()

L1.test <- subset(x=as.data.frame(layer1[1][[1]]),
       subset=(row.names(layer1[1][[1]]) %in% test.common))

L2.test <- subset(x=as.data.frame(layer2[1][[1]]),
             subset=(row.names(layer2[1][[1]]) %in% test.common))

# Então, gera-se a matriz final
profiles.test[[1]] = cbind(L1.test, L2.test, "Tag"=NA)

# Por fim, adiciona-se a nomenclatura com base nas colunas presentes
tag <- 0
if ("LFC.TS3.TN3" %in% colnames(P1.test)) {tag <- tag + 10000}
if ("LFC.TS6.TN6" %in% colnames(P1.test)) {tag <- tag + 100}
if ("LFC.TS9.TN9" %in% colnames(P1.test)) {tag <- tag + 1}
if ("L2FC.S.60.30" %in% colnames(P1.test)) {tag <- tag + 2000}
if ("L2FC.N.60.30" %in% colnames(P1.test)) {tag <- tag + 3000}
if ("L2FC.S.90.60" %in% colnames(P1.test)) {tag <- tag + 20}
if ("L2FC.N.90.60" %in% colnames(P1.test)) {tag <- tag + 30}

this.lt <- length(profiles.test[[1]])
profiles.test[[1]][this.lt] <- tag

## Por fim, confere se está tudo certo:
profiles.test[[1]]


##### Laço concatenando os 91 perfis ######
# 1. Primeiro, descobre-se o que eles tem em comum
# 2. Então, usando este nome, inicializa-se na lista de perfis
# 3. Cola-se a matriz final
# 4. Adiciona-se a tag

# Inicializando-se a lista de perfis
profiles <- list()
# Inicializa o número de perfis
n_prof <- 1

# Para cada camada externa,
for (L1.layer in layer1[1:length(layer1)]) {
  # Percorrendo cada camada interna,
  for (L2.layer in layer2[1:length(layer2)]) {
    
    # Auto-test 1
    # print(length(L1.layer))
    # print(length(L2.layer))
    
    # Pegue o nome de todas as isoformas em comum entre elas
    profile.names <- subset(
      x=row.names(L1.layer),
      row.names(L1.layer) %in% row.names(L2.layer)
      )
    
      ## ERROR CONTROL
      # Só continua rodando se houver algo em comum
      if (length(profile.names) != 0) {
    
    # Com os nomes em mãos, volte nas matrizes de DE para pegar as informações das isoformas
      ## Primeiro da camada 01 (S30/N30, S60/N60, S90/N90)
      info.L1 <- subset(
        x=as.data.frame(L1.layer),
        subset=(row.names(L1.layer) %in% profile.names)
      )
    
      ## Agora, a camada 02 (TC60/30, TC90/60, para S e N)
      info.L2 <- subset(
        x=as.data.frame(L2.layer),
        subset=(row.names(L2.layer) %in% profile.names)
      )
      
      # ERROR CONTROL
      # Verifica se ambos são do mesmo tamanho
      if (length(row.names(info.L2)) != length(row.names(info.L1))) {
      print("Something wrong happened!")
      }
      
      # Combina as duas camadas
      profiles[[n_prof]] = cbind(info.L1, info.L2, "Tag"=NA)
      
      # Por fim, adiciona-se a nomenclatura com base nas colunas presentes
      tag <- 0
      if ("LFC.TS3.TN3" %in% colnames(profiles[[n_prof]])) {tag <- tag + 10000}
      if ("LFC.TS6.TN6" %in% colnames(profiles[[n_prof]])) {tag <- tag + 100}
      if ("LFC.TS9.TN9" %in% colnames(profiles[[n_prof]])) {tag <- tag + 1}
      if ("L2FC.S.60.30" %in% colnames(profiles[[n_prof]])) {tag <- tag + 2000}
      if ("L2FC.N.60.30" %in% colnames(profiles[[n_prof]])) {tag <- tag + 3000}
      if ("L2FC.S.90.60" %in% colnames(profiles[[n_prof]])) {tag <- tag + 20}
      if ("L2FC.N.90.60" %in% colnames(profiles[[n_prof]])) {tag <- tag + 30}
      
      ## Adiciona-se os zeros à tag
      tag <- str_pad(tag, 5, pad = "0")
      
      ## Cola a tag na nova matriz criada
      profile.size <- length(profiles[[n_prof]])
      profiles[[n_prof]][profile.size] <- tag
    
      ## Soma 1 ao número de perfis
      n_prof <- n_prof + 1
    }
  }
}

## Descomente para zerar tudo, e começar novamente
# rm(profiles)
# rm(profile.size)

# De 91 perfis preditos, 74 existem:
length(profiles)

###### Passo 4. Adicionando se as isoformas estão ou não fragmentadas ######

# Com os perfis em mãos, agora adicionaremos a informação do TransRate,
# se as isoformas estão bem montadas (Complete) ou fragmentadas (Fragmented).

###### Passo 5. Adicionando anotações a cada isoforma ######

###### Passo 6. Exportando cada perfil como tabela com o nome de sua tag ######

##### WIP: Criando os perfis da planta suscetível ######
# Deixando para um momento futuro, quando essa informação for conveniente.

