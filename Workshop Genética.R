library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(data.table)
library(Matrix)
library(dplyr)
library(readxl)
library(wrapr)

setwd("~/Library/CloudStorage/GoogleDrive-caio.mussatto@unesp.br/Meu Drive/Projeto_de_Pesquisa /DATA")

linhagens <- readRDS('linhagens.rds')

# Lendo os CIFS

lista_CIFS <- read.table('CIFS_filtrada.txt')

# Lung

Lung <- subset(linhagens, subset = source == 'LUNG')

Lung <- RunUMAP(Lung, dims = 1:15, umap.method = 'umap-learn', metric = 'correlation')

Dimplot <- DimPlot(Lung, reduction = 'umap', group.by = 'sample')

ggsave('Imagens/Lung_Dimplot.tiff', plot = Dimplot, width = 10, height = 10)


dot <- DotPlot(Lung, features = lista_CIFS)
head(dot)

df <- dot$data

df <- df[, 5:4]

id <- df['id']

id <- distinct(id)

avg_lung <- data.frame(matrix(NA, nrow = 102, ncol = 5))

rownames(avg_lung) <- lista_CIFS$V1

colnames(avg_lung) <- id$id

avg_lung$NCIH2347_LUNG <- filter(df, id == 'NCIH2347_LUNG')
avg_lung$RERFLCAD1_LUNG <- filter(df, id == 'RERFLCAD1_LUNG')
avg_lung$RERFLCKJ_LUNG <- filter(df, id == 'RERFLCKJ_LUNG')
avg_lung$NCIH2087_LUNG <- filter(df, id == 'NCIH2087_LUNG')
avg_lung$NCIH358_LUNG <- filter(df, id == 'NCIH358_LUNG')

write.csv(avg_lung, 'Morpheus_avg_lung.csv')
lista_CIFS_lung <- read_csv("Workshop genetica /Morpheus_avg_lung.csv")

i = 3

while (i <= 11){
  lista_CIFS_lung <- lista_CIFS_lung[, -i]
  i <- i + 1
}

lista_CIFS_lung <- na.omit(lista_CIFS_lung)

VlnPlot1 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[1:12])
VlnPlot2 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[13:24])
VlnPlot3 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[25:36])
VlnPlot4 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[37:48])
VlnPlot5 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[49:60])
VlnPlot6 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[61:72])
VlnPlot7 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[73:81])
VlnPlot8 <- VlnPlot(Lung, features = lista_CIFS_lung$...1[82:86])

ggsave('Imagens/VLnPLot1.tiff', plot = VlnPlot1, width = 10, height = 10)
ggsave('Imagens/VlnPlot2.tiff', plot = VlnPlot2, width = 10, height = 10)
ggsave('Imagens/VlnPlot3.tiff', plot = VlnPlot3, width = 10, height = 10)
ggsave('Imagens/VlnPlot4.tiff', plot = VlnPlot4, width = 10, height = 10)
ggsave('Imagens/VlnPlot5.tiff', plot = VlnPlot5, width = 10, height = 10)
ggsave('Imagens/VlnPlot6.tiff', plot = VlnPlot6, width = 10, height = 10)
ggsave('Imagens/VlnPlot7.tiff', plot = VlnPlot7, width = 10, height = 10)
ggsave('Imagens/VlnPlot8.tiff', plot = VlnPlot8, width = 10, height = 10)


FeaturePlot1 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[1:12])
FeaturePlot2 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[13:24])
FeaturePlot3 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[25:36])
FeaturePlot4 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[37:48])
FeaturePlot5 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[49:60])
FeaturePlot6 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[61:72])
FeaturePlot7 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[73:81])
FeaturePlot8 <- FeaturePlot(Lung, features = lista_CIFS_lung$...1[82:86])


ggsave('Imagens/FeaturePlot1.tiff', plot = FeaturePlot1, width = 10, height = 10)
ggsave('Imagens/FeaturePlot2.tiff', plot = FeaturePlot2, width = 10, height = 10)
ggsave('Imagens/FeaturePlot3.tiff', plot = FeaturePlot3, width = 10, height = 10)
ggsave('Imagens/FeaturePlot4.tiff', plot = FeaturePlot4, width = 10, height = 10)
ggsave('Imagens/FeaturePlot5.tiff', plot = FeaturePlot5, width = 10, height = 10)
ggsave('Imagens/FeaturePlot6.tiff', plot = FeaturePlot6, width = 10, height = 10)
ggsave('Imagens/FeaturePlot7.tiff', plot = FeaturePlot7, width = 10, height = 10)
ggsave('Imagens/FeaturePlot8.tiff', plot = FeaturePlot8, width = 10, height = 10)

# Pancreas 

Pancreas <- subset(linhagens, subset = source == 'PANCREAS')

Pancreas <- RunUMAP(Pancreas, dims = 1:15, umap.method = 'umap-learn', metric = 'correlation')

Dimplot <- DimPlot(Pancreas, reduction = 'umap', group.by = 'sample')

ggsave('Imagens/Pancreas_Dimplot.tiff', plot = Dimplot, width = 10, height = 10)


dot <- DotPlot(Pancreas, features = lista_CIFS)
head(dot)

df <- dot$data

df <- df[, 5:4]

id <- df['id']

id <- distinct(id)

avg_Pancreas <- data.frame(matrix(NA, nrow = 102, ncol = 5))

rownames(avg_Pancreas) <- lista_CIFS$V1

colnames(avg_Pancreas) <- id$id

avg_Pancreas$HUPT3_PANCREAS <- filter(df, id == 'HUPT3_PANCREAS')
avg_Pancreas$HUPT4_PANCREAS <- filter(df, id == 'HUPT4_PANCREAS')
avg_Pancreas$PANC0203_PANCREAS <- filter(df, id == 'PANC0203_PANCREAS')
avg_Pancreas$ASPC1_PANCREAS <- filter(df, id == 'ASPC1_PANCREAS')
avg_Pancreas$SU8686_PANCREAS <- filter(df, id == 'SU8686_PANCREAS')

write.csv(avg_Pancreas, 'Morpheus_avg_Pancreas.csv')
lista_CIFS_Pancreas <- read_csv("Workshop genetica /Morpheus_avg_Pancreas.csv")

i = 3

while (i <= 11){
  lista_CIFS_Pancreas <- lista_CIFS_Pancreas[, -i]
  i <- i + 1
}

lista_CIFS_Pancreas <- na.omit(lista_CIFS_Pancreas)

VlnPlot1 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[1:12])
VlnPlot2 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[13:24])
VlnPlot3 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[25:36])
VlnPlot4 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[37:48])
VlnPlot5 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[49:60])
VlnPlot6 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[61:72])
VlnPlot7 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[73:81])
VlnPlot8 <- VlnPlot(Pancreas, features = lista_CIFS_Pancreas$...1[82:87])

ggsave('Imagens/VLnPLot1.tiff', plot = VlnPlot1, width = 10, height = 10)
ggsave('Imagens/VlnPlot2.tiff', plot = VlnPlot2, width = 10, height = 10)
ggsave('Imagens/VlnPlot3.tiff', plot = VlnPlot3, width = 10, height = 10)
ggsave('Imagens/VlnPlot4.tiff', plot = VlnPlot4, width = 10, height = 10)
ggsave('Imagens/VlnPlot5.tiff', plot = VlnPlot5, width = 10, height = 10)
ggsave('Imagens/VlnPlot6.tiff', plot = VlnPlot6, width = 10, height = 10)
ggsave('Imagens/VlnPlot7.tiff', plot = VlnPlot7, width = 10, height = 10)
ggsave('Imagens/VlnPlot8.tiff', plot = VlnPlot8, width = 10, height = 10)


FeaturePlot1 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[1:12])
FeaturePlot2 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[13:24])
FeaturePlot3 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[25:36])
FeaturePlot4 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[37:48])
FeaturePlot5 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[49:60])
FeaturePlot6 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[61:72])
FeaturePlot7 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[73:81])
FeaturePlot8 <- FeaturePlot(Pancreas, features = lista_CIFS_Pancreas$...1[82:87])


ggsave('Imagens/FeaturePlot1.tiff', plot = FeaturePlot1, width = 10, height = 10)
ggsave('Imagens/FeaturePlot2.tiff', plot = FeaturePlot2, width = 10, height = 10)
ggsave('Imagens/FeaturePlot3.tiff', plot = FeaturePlot3, width = 10, height = 10)
ggsave('Imagens/FeaturePlot4.tiff', plot = FeaturePlot4, width = 10, height = 10)
ggsave('Imagens/FeaturePlot5.tiff', plot = FeaturePlot5, width = 10, height = 10)
ggsave('Imagens/FeaturePlot6.tiff', plot = FeaturePlot6, width = 10, height = 10)
ggsave('Imagens/FeaturePlot7.tiff', plot = FeaturePlot7, width = 10, height = 10)
ggsave('Imagens/FeaturePlot8.tiff', plot = FeaturePlot8, width = 10, height = 10)