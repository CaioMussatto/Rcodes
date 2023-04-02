library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(data.table)
library(Matrix)
library(dplyr)
library(readxl)
library(wrapr)

getwd()


setwd("~/Library/CloudStorage/GoogleDrive-caio.mussatto@unesp.br/Meu Drive/Projeto_de_Pesquisa /DATA")

# Lendo os arquivos

pre_seurat <- readMM('Exp_data_TPM.mtx')

genes <- read_table('Genes.txt', col_names = FALSE)
gene.name <- genes[["X1"]]
cell_ids <- read.csv('Cells.csv')
cell.name <- cell_ids[["cell_name"]]


# Adicionando um filtro nas celulas 

metadata2 <- read_csv('Metadata.csv')

cells <- cell_ids[which(cell_ids$sample %in% metadata2$Cell_line),]

cells[['Filtro']] <- 2

cells <- left_join(cell_ids, cells)

cells[is.na(cells)] <- 1

cell_filtro <- cells[['Filtro']]

cell_filtro <- as.data.frame(cell_filtro)

rownames(cell_filtro) <- cells[['cell_name']]

2 %in% cells$Filtro

1 %in% cells$Filtro

# Criando o objeto Seurat

rownames(pre_seurat) <- gene.name
colnames(pre_seurat) <- cell.name


linhagens <- CreateSeuratObject(counts = pre_seurat, project = 'Linhagens',  min.cells = 3, min.features = 20)

sample <- cell_ids[['sample']]

source <- cell_ids[['source']]

linhagens$sample <- sample

linhagens$source <- source

# Realizando a filtragem 

linhagens <- AddMetaData(linhagens, cell_filtro)

linhagens <- subset(linhagens, subset = cell_filtro > 1 )

2 %in% linhagens$cell_filtro

1 %in% linhagens$cell_filtro

linhagens$cell_filtro <- NULL

# Lendo o arquivo RDS 

linhagens <- readRDS('Seurat/linhagens.rds')

View(linhagens)

# Normalization 
linhagens <- NormalizeData(linhagens)

# Feature select 

linhagens <- FindVariableFeatures(object = linhagens, selection.method = "vst", nfeatures = 2000)

top <- head(VariableFeatures(linhagens), 10)

plot1 <- VariableFeaturePlot(linhagens)

plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)

plot3 = plot1 + plot2

ggsave('Seurat/Imagens/FeatureSelect.tiff', plot = plot3, width = 10, height = 10)

#Scaling the data

todos_genes <- rownames(linhagens)

linhagens <- ScaleData(linhagens, features = todos_genes)

#Linear Reduction

linhagens <- RunPCA(linhagens, features = VariableFeatures(object = linhagens))

print(linhagens[["pca"]], dims = 1:5, nfeatures = 5)

# Ou

VizDim <- VizDimLoadings(linhagens, dims = 1:2, reduction = "pca")

ggsave('Seurat/Imagens/VizDim.tiff', plot = VizDim, width = 10, height = 10)


# Ou

DimPlot(linhagens, reduction = "pca")

DimHeatmap(linhagens, dims = 1, cells = 500, balanced = TRUE)

a <- DimHeatmap(linhagens, dims = 1:10, cells = 500, balanced = TRUE)

ggsave('/Imagens/DimHeatmap.tiff', plot = DimHeatmap, width = 10, height = 10)

linhagens <- readRDS('Seurat/linhagens.rds')

linhagens <- JackStraw(linhagens, num.replicate = 100)
linhagens <- ScoreJackStraw(linhagens, dims = 1:15)

JackStrawPlot(linhagens, dims = 1:15)
ElbowPlot <- ElbowPlot(linhagens)

ggsave('/Imagens/Elbowplot.tiff', plot = ElbowPlot, width = 10, height = 10)


# Clustering

linhagens <- FindNeighbors(linhagens, dims = 1:15)
linhagens <- FindClusters(linhagens, resolution = 0.5)


linhagens <- readRDS('linhagens_scaling.rds')

head(Idents(linhagens), 5)
                         
linhagens <- RunUMAP(linhagens, dims = 1:15, umap.method = 'umap-learn', metric = 'correlation')

Dimplot <- DimPlot(linhagens, reduction = 'umap', group.by = 'source')

ggsave('Imagens/Dimplot.tiff', plot = Dimplot, width = 10, height = 10)


linhagens <- readRDS('linhagens.rds')

Idents(linhagens) <- linhagens$sample

linhagens <- readRDS('linhagens.rds')

#lista_CIFS <- read.table('CIFS.txt')

lista_CIFS <- read.table('CIFS_filtrada.txt')

#teste <- anti_join(lista_CIFS, lista_CIFS_filtrada, by = 'V1')

#write_csv(teste, 'Genes_excluidos.csv')

# Dotplot de todos as linhagens

dot <- DotPlot(linhagens, features = lista_CIFS)

head(dot)

df <- dot$data


genes <- as.data.frame(df$features.plot)

genes <- distinct(genes)

colnames(genes) <- 'V1'

lista_CIFS <- inner_join(lista_CIFS, genes)


df <- read.csv('Linhagens_CIFS.csv')

df <- df[, 5:4]

id <- df['id']

id <- distinct(id)

avg <- data.frame(matrix(NA, nrow = 102, ncol = 84))

rownames(avg) <- lista_CIFS$V1


colnames(avg) <- id$id



avg$SW579_THYROID <- filter(df, id == 'SW579_THYROID')
avg$HEC251_ENDOMETRIUM <- filter(df, id == 'HEC251_ENDOMETRIUM')
avg$MFE319_ENDOMETRIUM <- filter(df, id == 'MFE319_ENDOMETRIUM')
avg$COLO741_SKIN <- filter(df, id == 'COLO741_SKIN')
avg$JHH7_LIVER <- filter(df, id == 'JHH7_LIVER')
avg$KNS42_CENTRAL_NERVOUS_SYSTEM <- filter(df, id == 'KNS42_CENTRAL_NERVOUS_SYSTEM')
avg$HT1197_URINARY_TRACT <- filter(df, id == 'HT1197_URINARY_TRACT')
avg$HT1376_URINARY_TRACT <- filter(df, id == 'HT1376_URINARY_TRACT')
avg$SNU308_BILIARY_TRACT <- filter(df, id == 'SNU308_BILIARY_TRACT')
avg$TM31_CENTRAL_NERVOUS_SYSTEM <- filter(df, id == 'TM31_CENTRAL_NERVOUS_SYSTEM')
avg$BT474_BREAST <- filter(df, id == 'BT474_BREAST')
avg$DKMG_CENTRAL_NERVOUS_SYSTEM <- filter(df, id == 'DKMG_CENTRAL_NERVOUS_SYSTEM')
avg$BT549_BREAST <- filter(df, id == 'BT549_BREAST')
avg$BICR6_UPPER_AERODIGESTIVE_TRACT <- filter(df, id == 'BICR6_UPPER_AERODIGESTIVE_TRACT')
avg$SH10TC_STOMACH <- filter(df, id == 'SH10TC_STOMACH')
avg$UMUC1_URINARY_TRACT <- filter(df, id == 'UMUC1_URINARY_TRACT')
avg$LS1034_LARGE_INTESTINE <- filter(df, id == 'LS1034_LARGE_INTESTINE')
avg$CCFSTTG1_CENTRAL_NERVOUS_SYSTEM <- filter(df, id == 'CCFSTTG1_CENTRAL_NERVOUS_SYSTEM')
avg$LNCAPCLONEFGC_PROSTATE <- filter(df, id == 'LNCAPCLONEFGC_PROSTATE')
avg$RCC10RGB_KIDNEY <- filter(df, id == 'RCC10RGB_KIDNEY')
avg$NCIH2347_LUNG <- filter(df, id == 'NCIH2347_LUNG')
avg$RERFLCAD1_LUNG <- filter(df, id == 'RERFLCAD1_LUNG')
avg$BICR31_UPPER_AERODIGESTIVE_TRACT <- filter(df, id == 'BICR31_UPPER_AERODIGESTIVE_TRACT')
avg$SKMEL3_SKIN <- filter(df, id == 'SKMEL3_SKIN')
avg$SNU1079_BILIARY_TRACT <- filter(df, id == 'SNU1079_BILIARY_TRACT')
avg$IGR1_SKIN <- filter(df, id == 'IGR1_SKIN')
avg$VMCUB1_URINARY_TRACT <- filter(df, id == 'VMCUB1_URINARY_TRACT')
avg$`42MGBA_CENTRAL_NERVOUS_SYSTEM` <- filter(df, id == '42MGBA_CENTRAL_NERVOUS_SYSTEM')
avg$LI7_LIVER <- filter(df, id == 'LI7_LIVER')
avg$FTC133_THYROID <- filter(df, id == 'FTC133_THYROID')
avg$SJSA1_BONE <- filter(df, id == 'SJSA1_BONE')
avg$HUPT3_PANCREAS <- filter(df, id == 'HUPT3_PANCREAS')
avg$BCPAP_THYROID <- filter(df, id == 'BCPAP_THYROID')
avg$PECAPJ49_UPPER_AERODIGESTIVE_TRACT <- filter(df, id == 'PECAPJ49_UPPER_AERODIGESTIVE_TRACT')
avg$`2313287_STOMACH` <- filter(df, id == '2313287_STOMACH')
avg$OSRC2_KIDNEY <- filter(df, id == 'OSRC2_KIDNEY')
avg$RERFLCKJ_LUNG <- filter(df, id == 'RERFLCKJ_LUNG')
avg$HUPT4_PANCREAS <- filter(df, id == 'HUPT4_PANCREAS')
avg$VMRCRCZ_KIDNEY <- filter(df, id == 'VMRCRCZ_KIDNEY')
avg$HUH6_LIVER <- filter(df, id == 'HUH6_LIVER')
avg$RVH421_SKIN <- filter(df, id == 'RVH421_SKIN')
avg$OVSAHO_OVARY <- filter(df, id == 'OVSAHO_OVARY')
avg$MFE280_ENDOMETRIUM <- filter(df, id == 'MFE280_ENDOMETRIUM')
avg$PC3_PROSTATE <- filter(df, id == 'PC3_PROSTATE')
avg$JHOS2_OVARY <- filter(df, id == 'JHOS2_OVARY')
avg$TE1_OESOPHAGUS <- filter(df, id == 'TE1_OESOPHAGUS')
avg$NCIH2087_LUNG <- filter(df, id == 'NCIH2087_LUNG')
avg$HEC59_ENDOMETRIUM <- filter(df, id == 'HEC59_ENDOMETRIUM')
avg$EFM192A_BREAST <- filter(df, id == 'EFM192A_BREAST')
avg$NCIH358_LUNG <- filter(df, id == 'NCIH358_LUNG')
avg$KYSE520_OESOPHAGUS <- filter(df, id == 'KYSE520_OESOPHAGUS')
avg$KMRC3_KIDNEY <- filter(df, id == 'KMRC3_KIDNEY')
avg$MKN7_STOMACH <- filter(df, id == 'MKN7_STOMACH')
avg$HCC1419_BREAST <- filter(df, id == 'HCC1419_BREAST')
avg$PANC0203_PANCREAS <- filter(df, id == 'PANC0203_PANCREAS')
avg$CL34_LARGE_INTESTINE <- filter(df, id == 'CL34_LARGE_INTESTINE')
avg$LS180_LARGE_INTESTINE <- filter(df, id == 'LS180_LARGE_INTESTINE')
avg$CAMA1_BREAST <- filter(df, id == 'CAMA1_BREAST')
avg$SNU449_LIVER <- filter(df, id == 'SNU449_LIVER')
avg$HUCCT1_BILIARY_TRACT <- filter(df, id == 'HUCCT1_BILIARY_TRACT')
avg$MKN45_STOMACH <- filter(df, id == 'MKN45_STOMACH')
avg$SNU1196_BILIARY_TRACT <- filter(df, id == 'SNU1196_BILIARY_TRACT')
avg$HEP3B217_LIVER <- filter(df, id == 'HEP3B217_LIVER')
avg$HOS_BONE <- filter(df, id == 'HOS_BONE')
avg$HEC151_ENDOMETRIUM <- filter(df, id == 'HEC151_ENDOMETRIUM')
avg$TE9_OESOPHAGUS <- filter(df, id == 'TE9_OESOPHAGUS')
avg$HCC56_LARGE_INTESTINE <- filter(df, id == 'HCC56_LARGE_INTESTINE')
avg$HT55_LARGE_INTESTINE <- filter(df, id == 'HT55_LARGE_INTESTINE')
avg$RT4_URINARY_TRACT <- filter(df, id == 'RT4_URINARY_TRACT')
avg$SKES1_BONE <- filter(df, id == 'SKES1_BONE')
avg$`8305C_THYROID` <- filter(df, id == '8305C_THYROID')
avg$TE10_OESOPHAGUS <- filter(df, id == 'TE10_OESOPHAGUS')
avg$ASPC1_PANCREAS <- filter(df, id == 'ASPC1_PANCREAS')
avg$BICR16_UPPER_AERODIGESTIVE_TRACT <- filter(df, id == 'BICR16_UPPER_AERODIGESTIVE_TRACT')
avg$OAW28_OVARY <- filter(df, id == 'OAW28_OVARY')
avg$SKMEL5_SKIN <- filter(df, id == 'SKMEL5_SKIN')
avg$OVCAR4_OVARY <- filter(df, id == 'OVCAR4_OVARY')
avg$ONCODG1_OVARY <- filter(df, id == 'ONCODG1_OVARY')
avg$HUH28_BILIARY_TRACT <- filter(df, id == 'HUH28_BILIARY_TRACT')
avg$IM95_STOMACH <- filter(df, id == 'IM95_STOMACH')
avg$TE14_OESOPHAGUS <- filter(df, id == 'TE14_OESOPHAGUS')
avg$CAKI2_KIDNEY <- filter(df, id == 'CAKI2_KIDNEY')
avg$SU8686_PANCREAS <- filter(df, id == 'SU8686_PANCREAS')
avg$SCC25_UPPER_AERODIGESTIVE_TRACT <- filter(df, id == 'SCC25_UPPER_AERODIGESTIVE_TRACT')



write.csv(avg, 'Morpheus_avg.csv')

avg <- read_csv('Morpheus_avg.csv')

i = 3

while (i <= 169){
  avg <- avg[, -i]
  i <- i + 1
}



write_csv(avg, 'Morpheus_avg_filtrado.csv')

avg <- read.csv('Morpheus_avg_filtrado.csv')

rownames(avg) <- avg[,1]

avg <- avg[-1,]

rownames(avg) <- avg$...1

avg <- avg[,-1]

linhagens_sample <- as_tibble(linhagens$sample)

linhagens_sample <- distinct(linhagens_sample)

write_csv(linhagens_sample, 'linhagens_sample.csv')

linhagens_sample <- read_csv('linhagens_sample.csv')


