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

linhagens_sample <- read_csv('linhagens_sample.csv')

SW579_THYROID.markers <- FindMarkers(linhagens, ident.1 = linhagens_sample$value[1], ident.2 = linhagens_sample$value[c(2,3,4)])

write.csv(SW579_THYROID.markers, 'Markers/SW579_THYROID.markers.csv')

FTC133_THYROID.markers <- FindMarkers(linhagens, ident.1 = linhagens_sample$value[2], ident.2 = linhagens_sample$value[c(1,3,4)])

write.csv(FTC133_THYROID.markers, 'Markers/FTC133_THYROID.markers.csv')

BCPAP_THYROID.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[3], ident.2 = linhagens_sample$value[c(1,2,4)])

write.csv(BCPAP_THYROID.markers, 'Markers/BCPAP_THYROID.markers.csv')

`8305C_THYROID.markers` <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[4], ident.2 = linhagens_sample$value[c(1,2,3)])

write.csv(`8305C_THYROID.markers`, 'Markers/8305C_THYROID.markers.csv')

HEC251_ENDOMETRIUM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[5], ident.2 = linhagens_sample$value[c(6,7,8,9)])

write.csv(HEC251_ENDOMETRIUM.markers, 'Markers/HEC251_ENDOMETRIUM.markers.csv')

MFE319_ENDOMETRIUM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[6], ident.2 = linhagens_sample$value[c(5,7,8,9)])

write.csv(MFE319_ENDOMETRIUM.markers, 'Markers/MFE319_ENDOMETRIUM.markers.csv')

MFE280_ENDOMETRIUM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[7], ident.2 = linhagens_sample$value[c(5,6,8,9)])

write.csv(MFE280_ENDOMETRIUM.markers, 'Markers/MFE280_ENDOMETRIUM.markers.csv')

HEC59_ENDOMETRIUM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[8], ident.2 = linhagens_sample$value[c(5,6,7,9)])

write.csv(HEC59_ENDOMETRIUM.markers, 'Markers/HEC59_ENDOMETRIUM.markers.csv')

HEC151_ENDOMETRIUM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[9], ident.2 = linhagens_sample$value[c(5,6,7,8)])

write.csv(HEC151_ENDOMETRIUM.markers, 'Markers/HEC151_ENDOMETRIUM.markers.csv')

COLO741_SKIN.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[10], ident.2 = linhagens_sample$value[c(11,12,13,14)])

write.csv(COLO741_SKIN.markers, 'Markers/COLO741_SKIN.markers.csv')

SKMEL3_SKIN.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[11], ident.2 = linhagens_sample$value[c(10,12,13,14)])

write.csv(SKMEL3_SKIN.markers, 'Markers/SKMEL3_SKIN.markers.csv')

IGR1_SKIN.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[12], ident.2 = linhagens_sample$value[c(10,11,13,14)])

write.csv(IGR1_SKIN.markers, 'Markers/IGR1_SKIN.markers.csv')

RVH421_SKIN.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[13], ident.2 = linhagens_sample$value[c(10,11,12,14)])

write.csv(RVH421_SKIN.markers, 'Markers/RVH421_SKIN.markers.csv')

SKMEL5_SKIN.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[14], ident.2 = linhagens_sample$value[c(10,11,12,13)])

write.csv(SKMEL5_SKIN.markers, 'Markers/SKMEL5_SKIN.markers.csv')

JHH7_LIVER.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[15], ident.2 = linhagens_sample$value[c(16,17,18,19)])

write.csv(JHH7_LIVER.markers, 'Markers/JHH7_LIVER.markers.csv')

LI7_LIVER.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[16], ident.2 = linhagens_sample$value[c(15,17,18,19)])

write.csv(LI7_LIVER.markers, 'Markers/LI7_LIVER.markers.csv')

HUH6_LIVER.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[17], ident.2 = linhagens_sample$value[c(15,16,18,19)])

write.csv(HUH6_LIVER.markers, 'Markers/HUH6_LIVER.markers.csv')

SNU449_LIVER.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[18], ident.2 = linhagens_sample$value[c(15,16,17,19)])

write.csv(SNU449_LIVER.markers, 'Markers/SNU449_LIVER.markers.csv')

HEP3B217_LIVER.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[19], ident.2 = linhagens_sample$value[c(15,16,17,18)])

write.csv(HEP3B217_LIVER.markers, 'Markers/HEP3B217_LIVER.markers.csv')

KNS42_CENTRAL_NERVOUS_SYSTEM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[20], ident.2 = linhagens_sample$value[c(21,22,23,24)])

write.csv(KNS42_CENTRAL_NERVOUS_SYSTEM.markers, 'Markers/KNS42_CENTRAL_NERVOUS_SYSTEM.markers.csv')

TM31_CENTRAL_NERVOUS_SYSTEM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[21], ident.2 = linhagens_sample$value[c(20,22,23,24)])

write.csv(TM31_CENTRAL_NERVOUS_SYSTEM.markers, 'Markers/TM31_CENTRAL_NERVOUS_SYSTEM.markers.csv')

DKMG_CENTRAL_NERVOUS_SYSTEM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[22], ident.2 = linhagens_sample$value[c(20,21,23,24)])

write.csv(DKMG_CENTRAL_NERVOUS_SYSTEM.markers, 'Markers/DKMG_CENTRAL_NERVOUS_SYSTEM.markers.csv')

CCFSTTG1_CENTRAL_NERVOUS_SYSTEM.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[23], ident.2 = linhagens_sample$value[c(20,21,22,24)])

write.csv(CCFSTTG1_CENTRAL_NERVOUS_SYSTEM.markers, 'Markers/CCFSTTG1_CENTRAL_NERVOUS_SYSTEM.markers.csv')

`42MGBA_CENTRAL_NERVOUS_SYSTEM.markers` <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[24], ident.2 = linhagens_sample$value[c(20,21,22,23)])

write.csv(`42MGBA_CENTRAL_NERVOUS_SYSTEM.markers`, 'Markers/42MGBA_CENTRAL_NERVOUS_SYSTEM.markers.csv')

HT1197_URINARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[25], ident.2 = linhagens_sample$value[c(26,27,28,29)])

write.csv(HT1197_URINARY_TRACT.markers, 'Markers/HT1197_URINARY_TRACT.markers.csv')

HT1376_URINARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[26], ident.2 = linhagens_sample$value[c(25,27, 28, 29)])

write.csv(HT1376_URINARY_TRACT.markers, 'Markers/HT1376_URINARY_TRACT.markers.csv')

UMUC1_URINARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[27], ident.2 = linhagens_sample$value[c(25,26,28, 29)])

write.csv(UMUC1_URINARY_TRACT.markers, 'Markers/UMUC1_URINARY_TRACT.markers.csv')

VMCUB1_URINARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[28], ident.2 = linhagens_sample$value[c(25,26, 27, 29)])

write.csv(VMCUB1_URINARY_TRACT.markers, 'Markers/VMCUB1_URINARY_TRACT.markers.csv')

RT4_URINARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[29], ident.2 = linhagens_sample$value[c(25,26,28, 27)])

write.csv(RT4_URINARY_TRACT.markers, 'Markers/RT4_URINARY_TRACT.markers.csv')

SNU308_BILIARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[30], ident.2 = linhagens_sample$value[c(31,32,33, 34)])

write.csv(SNU308_BILIARY_TRACT.markers, 'Markers/SNU308_BILIARY_TRACT.markers.csv')

SNU1079_BILIARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[31], ident.2 = linhagens_sample$value[c(30,32,33, 34)])

write.csv(SNU1079_BILIARY_TRACT.markers, 'Markers/SNU1079_BILIARY_TRACT.markers.csv')

HUCCT1_BILIARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[32], ident.2 = linhagens_sample$value[c(30,31,33,34)])

write.csv(HUCCT1_BILIARY_TRACT.markers, 'Markers/HUCCT1_BILIARY_TRACT.markers.csv')

SNU1196_BILIARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[33], ident.2 = linhagens_sample$value[c(30,31,32, 34)])

write.csv(SNU1196_BILIARY_TRACT.markers, 'Markers/SNU1196_BILIARY_TRACT.markers.csv')

HUH28_BILIARY_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[34], ident.2 = linhagens_sample$value[c(30,31,32, 33)])

write.csv(HUH28_BILIARY_TRACT.markers, 'Markers/HUH28_BILIARY_TRACT.markers.csv')

BT474_BREAST.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[35], ident.2 = linhagens_sample$value[c(36,37,38, 39)])

write.csv(BT474_BREAST.markers, 'Markers/BT474_BREAST.markers.csv')

BT549_BREAST.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[36], ident.2 = linhagens_sample$value[c(35,37,38, 39)])

write.csv(BT549_BREAST.markers, 'Markers/BT549_BREAST.markers.csv')

EFM192A_BREAST.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[37], ident.2 = linhagens_sample$value[c(35,36,38,39)])

write.csv(EFM192A_BREAST.markers, 'Markers/EFM192A_BREAST.markers.csv')

HCC1419_BREAST.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[38], ident.2 = linhagens_sample$value[c(35,36,37,39)])

write.csv(HCC1419_BREAST.markers, 'Markers/HCC1419_BREAST.markers.csv')

CAMA1_BREAST.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[39], ident.2 = linhagens_sample$value[c(35,36,37, 38)])

write.csv(CAMA1_BREAST.markers, 'Markers/CAMA1_BREAST.markers.csv')

BICR6_UPPER_AERODIGESTIVE_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[40], ident.2 = linhagens_sample$value[c(41,42,43,44)])

write.csv(BICR6_UPPER_AERODIGESTIVE_TRACT.markers, 'Markers/BICR6_UPPER_AERODIGESTIVE_TRACT.markers.csv')

BICR31_UPPER_AERODIGESTIVE_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[41], ident.2 = linhagens_sample$value[c(40,42,43,44)])

write.csv(BICR31_UPPER_AERODIGESTIVE_TRACT.markers, 'Markers/BICR31_UPPER_AERODIGESTIVE_TRACT.markers.csv')

PECAPJ49_UPPER_AERODIGESTIVE_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[42], ident.2 = linhagens_sample$value[c(40,41,43,44)])

write.csv(PECAPJ49_UPPER_AERODIGESTIVE_TRACT.markers, 'Markers/PECAPJ49_UPPER_AERODIGESTIVE_TRACT.markers.csv')

BICR16_UPPER_AERODIGESTIVE_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[43], ident.2 = linhagens_sample$value[c(40,41,42, 44)])

write.csv(BICR16_UPPER_AERODIGESTIVE_TRACT.markers, 'Markers/BICR16_UPPER_AERODIGESTIVE_TRACT.markers.csv')

SCC25_UPPER_AERODIGESTIVE_TRACT.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[44], ident.2 = linhagens_sample$value[c(40,41,42, 43)])

write.csv(SCC25_UPPER_AERODIGESTIVE_TRACT.markers, 'Markers/SCC25_UPPER_AERODIGESTIVE_TRACT.markers.markers.csv')

SH10TC_STOMACH.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[45], ident.2 = linhagens_sample$value[c(46,47,48, 49)])

write.csv(SH10TC_STOMACH.markers, 'Markers/SH10TC_STOMACH.markers.csv')

`2313287_STOMACH.markers` <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[46], ident.2 = linhagens_sample$value[c(45,47,48, 49)])

write.csv(`2313287_STOMACH.markers`, 'Markers/2313287_STOMACH.markers.csv')

MKN7_STOMACH.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[47], ident.2 = linhagens_sample$value[c(45,46,48, 49)])

write.csv(MKN7_STOMACH.markers, 'Markers/MKN7_STOMACH.markers.csv')

MKN45_STOMACH.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[48], ident.2 = linhagens_sample$value[c(45,46,47, 49)])

write.csv(MKN45_STOMACH.markers, 'Markers/MKN45_STOMACH.markers.csv')

IM95_STOMACH.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[49], ident.2 = linhagens_sample$value[c(45,46,47,48)])

write.csv(IM95_STOMACH.markers, 'Markers/IM95_STOMACH.markers.csv')

LS1034_LARGE_INTESTINE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[50], ident.2 = linhagens_sample$value[c(51,52,53, 54)])

write.csv(LS1034_LARGE_INTESTINE.markers, 'Markers/LS1034_LARGE_INTESTINE.markers.csv')

CL34_LARGE_INTESTINE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[51], ident.2 = linhagens_sample$value[c(50,52, 53, 54)])

write.csv(CL34_LARGE_INTESTINE.markers, 'Markers/CL34_LARGE_INTESTINE.markers.csv')

LS180_LARGE_INTESTINE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[52], ident.2 = linhagens_sample$value[c(50,51,53, 54)])

write.csv(LS180_LARGE_INTESTINE.markers, 'Markers/LS180_LARGE_INTESTINE.markers.csv')

HCC56_LARGE_INTESTINE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[53], ident.2 = linhagens_sample$value[c(50,51,53, 54)])

write.csv(HCC56_LARGE_INTESTINE.markers, 'Markers/HCC56_LARGE_INTESTINE.markers.csv')

HT55_LARGE_INTESTINE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[54], ident.2 = linhagens_sample$value[c(50,51,53, 54)])

write.csv(HT55_LARGE_INTESTINE.markers, 'Markers/HT55_LARGE_INTESTINE.markers.csv')

LNCAPCLONEFGC_PROSTATE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[55], ident.2 = linhagens_sample$value[c(56)])

write.csv(LNCAPCLONEFGC_PROSTATE.markers, 'Markers/LNCAPCLONEFGC_PROSTATE.markers.csv')

PC3_PROSTATE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[56], ident.2 = linhagens_sample$value[c(55)])

write.csv(PC3_PROSTATE.markers, 'Markers/PC3_PROSTATE.markers.csv')

RCC10RGB_KIDNEY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[57], ident.2 = linhagens_sample$value[c(58,59,60, 61)])

write.csv(RCC10RGB_KIDNEY.markers, 'Markers/RCC10RGB_KIDNEY.markers.csv')

OSRC2_KIDNEY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[58], ident.2 = linhagens_sample$value[c(57,59,60, 61)])

write.csv(OSRC2_KIDNEY.markers, 'Markers/OSRC2_KIDNEY.markers.csv')

VMRCRCZ_KIDNEY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[59], ident.2 = linhagens_sample$value[c(57,58,60, 61)])

write.csv(VMRCRCZ_KIDNEY.markers, 'Markers/VMRCRCZ_KIDNEY.markers.csv')

KMRC3_KIDNEY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[60], ident.2 = linhagens_sample$value[c(57,58,59, 61)])

write.csv(KMRC3_KIDNEY.markers, 'Markers/KMRC3_KIDNEY.markers.csv')

CAKI2_KIDNEY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[61], ident.2 = linhagens_sample$value[c(57,58,59, 60)])

write.csv(CAKI2_KIDNEY.markers, 'Markers/CAKI2_KIDNEY.markers.csv')

NCIH2347_LUNG.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[62], ident.2 = linhagens_sample$value[c(63,64,65, 66)])

write.csv(NCIH2347_LUNG.markers, 'Markers/NCIH2347_LUNG.markers.csv')

RERFLCAD1_LUNG.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[63], ident.2 = linhagens_sample$value[c(62,64,65, 66)])

write.csv(RERFLCAD1_LUNG.markers, 'Markers/RERFLCAD1_LUNG.markers.csv')

RERFLCKJ_LUNG.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[64], ident.2 = linhagens_sample$value[c(62,63,65, 66)])

write.csv(RERFLCKJ_LUNG.markers, 'Markers/RERFLCKJ_LUNG.markers.csv')

NCIH2087_LUNG.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[65], ident.2 = linhagens_sample$value[c(62,63,64, 66)])

write.csv(NCIH2087_LUNG.markers, 'Markers/NCIH2087_LUNG.markers.csv')

NCIH358_LUNG.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[66], ident.2 = linhagens_sample$value[c(62,63,64, 65)])

write.csv(NCIH358_LUNG.markers, 'Markers/NCIH358_LUNG.markers.csv')

SJSA1_BONE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[67], ident.2 = linhagens_sample$value[c(68,69)])

write.csv(SJSA1_BONE.markers, 'Markers/SJSA1_BONE.markers.csv')

HOS_BONE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[68], ident.2 = linhagens_sample$value[c(67,69)])

write.csv(HOS_BONE.markers, 'Markers/HOS_BONE.markers.csv')

SKES1_BONE.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[69], ident.2 = linhagens_sample$value[c(67,68)])

write.csv(SKES1_BONE.markers, 'Markers/SKES1_BONE.markers.csv')

HUPT3_PANCREAS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[70], ident.2 = linhagens_sample$value[c(71,72,73, 74)])

write.csv(HUPT3_PANCREAS.markers, 'Markers/HUPT3_PANCREAS.markers.csv')

HUPT4_PANCREAS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[71], ident.2 = linhagens_sample$value[c(70,72,73, 74)])

write.csv(HUPT4_PANCREAS.markers, 'Markers/HUPT4_PANCREAS.markers.csv')

PANC0203_PANCREAS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[72], ident.2 = linhagens_sample$value[c(70,71,73, 74)])

write.csv(PANC0203_PANCREAS.markers, 'Markers/PANC0203_PANCREAS.markers.csv')

ASPC1_PANCREAS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[73], ident.2 = linhagens_sample$value[c(70,71,72, 74)])

write.csv(ASPC1_PANCREAS.markers, 'Markers/ASPC1_PANCREAS.markers.csv')

SU8686_PANCREAS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[74], ident.2 = linhagens_sample$value[c(70,71,72,73)])

write.csv(SU8686_PANCREAS.markers, 'Markers/SU8686_PANCREAS.markers.csv')

OVSAHO_OVARY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[75], ident.2 = linhagens_sample$value[c(76,77,78, 79)])

write.csv(OVSAHO_OVARY.markers, 'Markers/OVSAHO_OVARY.markers.csv')

JHOS2_OVARY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[76], ident.2 = linhagens_sample$value[c(75,77,78, 79)])

write.csv(JHOS2_OVARY.markers, 'Markers/JHOS2_OVARY.markers.csv')

OAW28_OVARY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[77], ident.2 = linhagens_sample$value[c(75,76,78, 79)])

write.csv(OAW28_OVARY.markers, 'Markers/OAW28_OVARY.markers.csv')

OVCAR4_OVARY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[78], ident.2 = linhagens_sample$value[c(75,76,77, 79)])

write.csv(OVCAR4_OVARY.markers, 'Markers/OVCAR4_OVARY.markers.csv')

ONCODG1_OVARY.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[79], ident.2 = linhagens_sample$value[c(75,76,77, 78)])

write.csv(ONCODG1_OVARY.markers, 'Markers/ONCODG1_OVARY.markers.csv')

TE1_OESOPHAGUS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[80], ident.2 = linhagens_sample$value[c(81,82,83, 84)])

write.csv(TE1_OESOPHAGUS.markers, 'Markers/TE1_OESOPHAGUS.markers.csv')

KYSE520_OESOPHAGUS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[81], ident.2 = linhagens_sample$value[c(80,82,83, 84)])

write.csv(KYSE520_OESOPHAGUS.markers, 'Markers/KYSE520_OESOPHAGUS.markers.csv')

TE9_OESOPHAGUS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[82], ident.2 = linhagens_sample$value[c(80,81,83, 84)])

write.csv(TE9_OESOPHAGUS.markers, 'Markers/TE9_OESOPHAGUS.markers.csv')

TE10_OESOPHAGUS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[83], ident.2 = linhagens_sample$value[c(80,81,82, 84)])

write.csv(TE10_OESOPHAGUS.markers, 'Markers/TE10_OESOPHAGUS.markers.csv')

TE14_OESOPHAGUS.markers <-  FindMarkers(linhagens, ident.1 = linhagens_sample$value[84], ident.2 = linhagens_sample$value[c(80,81,82, 83)])

write.csv(TE14_OESOPHAGUS.markers, 'Markers/TE14_OESOPHAGUS.markers.csv')


levels(linhagens)

         