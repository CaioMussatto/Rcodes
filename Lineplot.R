library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(data.table)
library(Matrix)
library(dplyr)
library(readxl)
library(wrapr)
library(dplyr)
library(hrbrthemes)
library(viridis)

setwd("~/Library/CloudStorage/GoogleDrive-caio.mussatto@unesp.br/Meu Drive/Projeto_de_Pesquisa /DATA")

# Lendo os arquivos

linhagens <- readRDS('linhagens.rds')
lista_CIFS <- read.table('CIFS_filtrada.txt')

#linhagens_sample <- as_tibble(linhagens$sample)

#linhagens_sample <- distinct(linhagens_sample)

#write_csv(linhagens_sample, 'linhagens_sample.csv')

linhagens_sample <- read_csv('linhagens_sample.csv')

new_df <- data.frame(gene = c(1:102))

row.names(new_df) <- lista_CIFS$V1

count_matrix <- as.data.frame(GetAssayData(object = linhagens, slot = "counts")[paste(lista_CIFS$V1,sep=''),WhichCells(object = linhagens, ident = 'SW579_THYROID')])

colunas <- function(seurat, identa){
  identa <- as.data.frame(GetAssayData(object = seurat, slot = "counts")[paste(lista_CIFS$V1,sep=''),WhichCells(object = linhagens, ident = identa)])
  return(rowMeans(identa))
}

for (i in linhagens_sample$value){
  print(i)
  new_df[i] <- colunas(linhagens, i)
}

rownames(new_df) <- lista_CIFS$V1

new_df <- new_df[,-1]

write_csv(new_df, 'Media_dos_genes.csv')


new_df <- read_csv('Media_dos_genes.csv', spcs = TRUE)



new_df$gene = rownames(new_df)

corte <- new_df[,1:5]


for (i in colnames(new_df_loop)){
  print(max(new_df[i]))
  #b <- grafico2(new_df_loop, i)
  #lista[[i]] <- b
}


corte$gene <- rownames(new_df)

# Ver oms maiores valores 

new_df_loop <- new_df[,-1]



df_reshaped <- data.frame(x = corte$gene,                           
                          y = c(corte$SW579_THYROID, corte$FTC133_THYROID, corte$BCPAP_THYROID, corte$`8305C_THYROID`),
                          group = c(rep("SW579_THYROID", nrow(corte)),
                                    rep("FTC133_THYROID", nrow(corte)),
                                    rep("BCPAP_THYROID", nrow(corte)),
                                    rep('8305C_THYROID', nrow(corte))))


ggplot(df_reshaped, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')



endometrium <- data.frame(x = new_df$gene,
                          y = c(new_df$HEC251_ENDOMETRIUM, new_df$MFE319_ENDOMETRIUM, new_df$MFE280_ENDOMETRIUM, new_df$HEC59_ENDOMETRIUM, new_df$HEC151_ENDOMETRIUM),
                          group = c(rep('HEC251_ENDOMETRIUM', nrow(new_df)),
                                    rep('MFE319_ENDOMETRIUM', nrow(new_df)),
                                    rep('MFE280_ENDOMETRIUM', nrow(new_df)),
                                    rep('HEC59_ENDOMETRIUM', nrow(new_df)),
                                    rep('HEC151_ENDOMETRIUM', nrow(new_df))))

ggplot(endometrium, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')


skin <- data.frame(x = new_df$gene,
                   y = c(new_df$COLO741_SKIN, new_df$SKMEL3_SKIN, new_df$IGR1_SKIN, new_df$RVH421_SKIN, new_df$SKMEL5_SKIN),
                   group = c(rep('COLO741_SKIN', nrow(new_df)),
                             rep('SKMEL3_SKIN', nrow(new_df)),
                             rep('IGR1_SKIN', nrow(new_df)),
                             rep('RVH421_SKIN', nrow(new_df)),
                             rep('SKMEL5_SKIN', nrow(new_df))))


ggplot(endometrium, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')


liver <-  data.frame(x = new_df$gene,
                     y = c(new_df$JHH7_LIVER, new_df$LI7_LIVER, new_df$HUH6_LIVER, new_df$SNU449_LIVER, new_df$HEP3B217_LIVER),
                     group = c(rep('JHH7_LIVER', nrow(new_df)),
                               rep('LI7_LIVER', nrow(new_df)),
                               rep('HUH6_LIVER', nrow(new_df)),
                               rep('SNU449_LIVER', nrow(new_df)),
                               rep('HEP3B217_LIVER', nrow(new_df))))

ggplot(liver, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')

Central_nervous_system <-  data.frame(x = new_df$gene,
                     y = c(new_df$KNS42_CENTRAL_NERVOUS_SYSTEM, new_df$TM31_CENTRAL_NERVOUS_SYSTEM, new_df$DKMG_CENTRAL_NERVOUS_SYSTEM, new_df$CCFSTTG1_CENTRAL_NERVOUS_SYSTEM, new_df$`42MGBA_CENTRAL_NERVOUS_SYSTEM`),
                     group = c(rep('KNS42_CENTRAL_NERVOUS_SYSTEM', nrow(new_df)),
                               rep('TM31_CENTRAL_NERVOUS_SYSTEM', nrow(new_df)),
                               rep('DKMG_CENTRAL_NERVOUS_SYSTEM', nrow(new_df)),
                               rep('CCFSTTG1_CENTRAL_NERVOUS_SYSTEM', nrow(new_df)),
                               rep('42MGBA_CENTRAL_NERVOUS_SYSTEM', nrow(new_df))))



ggplot(Central_nervous_system, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')

urinary_tract <-  data.frame(x = new_df$gene,
                                      y = c(new_df$HT1197_URINARY_TRACT, new_df$HT1376_URINARY_TRACT, new_df$UMUC1_URINARY_TRACT, new_df$VMCUB1_URINARY_TRACT, new_df$RT4_URINARY_TRACT),
                                      group = c(rep('HT1197_URINARY_TRACT', nrow(new_df)),
                                                rep('HT1376_URINARY_TRACT', nrow(new_df)),
                                                rep('UMUC1_URINARY_TRACT', nrow(new_df)),
                                                rep('VMCUB1_URINARY_TRACT', nrow(new_df)),
                                                rep('RT4_URINARY_TRACT', nrow(new_df))))

ggplot(urinary_tract, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')


biliary_tract <-  data.frame(x = new_df$gene,
                             y = c(new_df$SNU308_BILIARY_TRACT, new_df$SNU1079_BILIARY_TRACT, new_df$HUCCT1_BILIARY_TRACT, new_df$SNU1196_BILIARY_TRACT, new_df$HUH28_BILIARY_TRACT),
                             group = c(rep('SNU308_BILIARY_TRACT', nrow(new_df)),
                                       rep('SNU1079_BILIARY_TRACT', nrow(new_df)),
                                       rep('HUCCT1_BILIARY_TRACT', nrow(new_df)),
                                       rep('SNU1196_BILIARY_TRACT', nrow(new_df)),
                                       rep('HUH28_BILIARY_TRACT', nrow(new_df))))

ggplot(biliary_tract, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')


breast <- data.frame(x = new_df$gene,
                     y = c(new_df$BT474_BREAST, new_df$BT549_BREAST, new_df$EFM192A_BREAST, new_df$HCC1419_BREAST, new_df$CAMA1_BREAST),
                     group = c(rep('BT474_BREAST', nrow(new_df)),
                               rep('BT549_BREAST', nrow(new_df)),
                               rep('EFM192A_BREAST', nrow(new_df)),
                               rep('HCC1419_BREAST', nrow(new_df)),
                               rep('CAMA1_BREAST', nrow(new_df))))

ggplot(breast, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('')


upper <- data.frame(x = new_df$gene,
                     y = c(new_df$BICR6_UPPER_AERODIGESTIVE_TRACT, new_df$BICR31_UPPER_AERODIGESTIVE_TRACT, new_df$PECAPJ49_UPPER_AERODIGESTIVE_TRACT, new_df$BICR16_UPPER_AERODIGESTIVE_TRACT, new_df$SCC25_UPPER_AERODIGESTIVE_TRACT),
                     group = c(rep('BICR6_UPPER_AERODIGESTIVE_TRACT', nrow(new_df)),
                               rep('BICR31_UPPER_AERODIGESTIVE_TRACT', nrow(new_df)),
                               rep('PECAPJ49_UPPER_AERODIGESTIVE_TRACT', nrow(new_df)),
                               rep('BICR16_UPPER_AERODIGESTIVE_TRACT', nrow(new_df)),
                               rep('SCC25_UPPER_AERODIGESTIVE_TRACT', nrow(new_df))))

ggplot(upper, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('UPPER_AERODIGESTIVE_TRACT') + xlab('')

stomach <- data.frame(x = new_df$gene,
                    y = c(new_df$SH10TC_STOMACH, new_df$`2313287_STOMACH`, new_df$MKN7_STOMACH, new_df$MKN45_STOMACH, new_df$IM95_STOMACH),
                    group = c(rep('SH10TC_STOMACH', nrow(new_df)),
                              rep('2313287_STOMACH', nrow(new_df)),
                              rep('MKN7_STOMACH', nrow(new_df)),
                              rep('MKN45_STOMACH', nrow(new_df)),
                              rep('IM95_STOMACH', nrow(new_df))))

ggplot(stomach, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Stomach') + xlab('')

large_intestine <- data.frame(x = new_df$gene,
                      y = c(new_df$LS1034_LARGE_INTESTINE, new_df$CL34_LARGE_INTESTINE, new_df$LS180_LARGE_INTESTINE, new_df$HCC56_LARGE_INTESTINE, new_df$HT55_LARGE_INTESTINE),
                      group = c(rep('LS1034_LARGE_INTESTINE', nrow(new_df)),
                                rep('CL34_LARGE_INTESTINE', nrow(new_df)),
                                rep('LS180_LARGE_INTESTINE', nrow(new_df)),
                                rep('HCC56_LARGE_INTESTINE', nrow(new_df)),
                                rep('HT55_LARGE_INTESTINE', nrow(new_df))))

ggplot(large_intestine, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Large Intestine') + xlab('')


Prostate <- data.frame(x = new_df$gene,
                              y = c(new_df$LNCAPCLONEFGC_PROSTATE, new_df$PC3_PROSTATE),
                              group = c(rep('LNCAPCLONEFGC_PROSTATE', nrow(new_df)),
                                        rep('PC3_PROSTATE', nrow(new_df))))

ggplot(Prostate, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Prostate') + xlab('')


Kidney <- data.frame(x = new_df$gene,
                     y = c(new_df$RCC10RGB_KIDNEY, new_df$OSRC2_KIDNEY, new_df$VMRCRCZ_KIDNEY, new_df$KMRC3_KIDNEY, new_df$CAKI2_KIDNEY),
                     group = c(rep('RCC10RGB_KIDNEY', nrow(new_df)),
                               rep('OSRC2_KIDNEY', nrow(new_df)),
                               rep('VMRCRCZ_KIDNEY', nrow(new_df)),
                               rep('KMRC3_KIDNEY', nrow(new_df)),
                               rep('CAKI2_KIDNEY', nrow(new_df))))

ggplot(Kidney, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Kidney') + xlab('')

Lung <- data.frame(x = new_df$gene,
                     y = c(new_df$NCIH2347_LUNG, new_df$RERFLCAD1_LUNG, new_df$RERFLCKJ_LUNG, new_df$NCIH2087_LUNG, new_df$NCIH358_LUNG),
                     group = c(rep('NCIH2347_LUNG', nrow(new_df)),
                               rep('RERFLCAD1_LUNG', nrow(new_df)),
                               rep('RERFLCKJ_LUNG', nrow(new_df)),
                               rep('NCIH2087_LUNG', nrow(new_df)),
                               rep('NCIH358_LUNG', nrow(new_df))))

ggplot(Lung, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Lung') + xlab('')


Bone <- data.frame(x = new_df$gene,
                   y = c(new_df$SJSA1_BONE, new_df$HOS_BONE, new_df$SKES1_BONE),
                   group = c(rep('SJSA1_BONE', nrow(new_df)),
                             rep('HOS_BONE', nrow(new_df)),
                             rep('SKES1_BONE', nrow(new_df))))

ggplot(Bone, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Bone') + xlab('')


Pancreas <- data.frame(x = new_df$gene,
                   y = c(new_df$HUPT3_PANCREAS, new_df$HUPT4_PANCREAS, new_df$PANC0203_PANCREAS, new_df$ASPC1_PANCREAS, new_df$SU8686_PANCREAS),
                   group = c(rep('HUPT3_PANCREAS', nrow(new_df)),
                             rep('HUPT4_PANCREAS', nrow(new_df)),
                             rep('PANC0203_PANCREAS', nrow(new_df)),
                             rep('ASPC1_PANCREAS', nrow(new_df)),
                             rep('SU8686_PANCREAS', nrow(new_df))))

ggplot(Pancreas, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Pancreas') + xlab('')

Ovary <- data.frame(x = new_df$gene,
                       y = c(new_df$OVSAHO_OVARY, new_df$JHOS2_OVARY, new_df$OAW28_OVARY, new_df$OVCAR4_OVARY, new_df$ONCODG1_OVARY),
                       group = c(rep('OVSAHO_OVARY', nrow(new_df)),
                                 rep('JHOS2_OVARY', nrow(new_df)),
                                 rep('OAW28_OVARY', nrow(new_df)),
                                 rep('OVCAR4_OVARY', nrow(new_df)),
                                 rep('ONCODG1_OVARY', nrow(new_df))))

ggplot(Ovary, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Ovary') + xlab('')

Oesophagus <- data.frame(x = new_df$gene,
                    y = c(new_df$TE1_OESOPHAGUS, new_df$KYSE520_OESOPHAGUS, new_df$TE9_OESOPHAGUS, new_df$TE10_OESOPHAGUS, new_df$TE14_OESOPHAGUS),
                    group = c(rep('TE1_OESOPHAGUS', nrow(new_df)),
                              rep('KYSE520_OESOPHAGUS', nrow(new_df)),
                              rep('TE9_OESOPHAGUS', nrow(new_df)),
                              rep('TE10_OESOPHAGUS', nrow(new_df)),
                              rep('TE14_OESOPHAGUS', nrow(new_df))))

ggplot(Oesophagus, aes(x, y, col = group, group = 1)) + geom_line() + facet_grid(group ~ .) + theme(text = element_text(size = 8),  axis.text.x = element_text(angle = 90, hjust = 1)) + ylim (0,3000) + ylab('') + ggtitle('Oesophagus') + xlab('')




for (i in colnames(new_df_loop)) {
  p <- grafico(new_df, i)
  #ggsave(filename = paste("plot_", i, ".png", sep = ""), plot = p)
}

for (i in colnames(new_df_loop)){
  print(i)
  }

p

a <- grafico(new_df_loop, colnames(new_df_loop))

a

for (i in colnames(new_df)){
  print(i)
  next
  print(i)
}