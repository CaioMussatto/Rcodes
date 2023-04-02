library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(data.table)
library(Matrix)
library(dplyr)
library(readxl)
library(wrapr)
library(skimr)
library(moderndive)
library(GGally)
library(smplot2)
library(corrplot)
library(RColorBrewer)
library(ggpubr)
library(ggridges)
library(ggplot2)
setwd("~/Library/CloudStorage/GoogleDrive-caio.mussatto@unesp.br/Meu Drive/Projeto_de_Pesquisa /DATA")



avg <- read.csv('Morpheus_avg_filtrado.csv', sep = ';')

avg <- na.omit(avg)


rownames(avg) <- avg$Gene

vi <- as.data.frame(read_xlsx('avg_exp_PAN.xlsx'))

colnames(vi)[1]  <- 'Gene'

vi <- vi %>% inner_join(avg, c("Gene"="Gene"))

rownames(vi) <- vi$Gene

vi <- vi %>% select(starts_with("Malignant"))

vi <- as.data.frame(vi)

THYROID <- avg %>% 
  select(starts_with("THYROID"))
colnames(THYROID) <- THYROID[1,] 
THYROID <- THYROID[-1,]

ENDOMETRIUM <- avg %>% 
  select(starts_with("ENDOMETRIUM"))
colnames(ENDOMETRIUM) <- ENDOMETRIUM[1,] 
ENDOMETRIUM <- ENDOMETRIUM[-1,]  

SKIN <- avg %>% 
  select(starts_with("SKIN"))
colnames(SKIN) <- SKIN[1,] 
SKIN <- SKIN[-1,]
SKIN <- as.data.frame(SKIN)
SKIN$malignant <- vi[,8]


LIVER <- avg %>% 
  select(starts_with("LIVER"))

colnames(LIVER) <- LIVER[1,] 
LIVER <- LIVER[-1,]

CENTRAL_NERVOUS_SYSTEM <- avg %>% 
  select(starts_with("BRAIN"))

colnames(CENTRAL_NERVOUS_SYSTEM) <- CENTRAL_NERVOUS_SYSTEM[1,] 
CENTRAL_NERVOUS_SYSTEM <- CENTRAL_NERVOUS_SYSTEM[-1,]

URINARY_TRACT <- avg %>% 
  select(starts_with("BLADDER"))

colnames(URINARY_TRACT) <- URINARY_TRACT[1,] 
URINARY_TRACT <- URINARY_TRACT[-1,]

BILIARY_TRACT <- avg %>% 
  select(starts_with("GALLBLADDER"))

colnames(BILIARY_TRACT) <- BILIARY_TRACT[1,] 
BILIARY_TRACT <- BILIARY_TRACT[-1,]

BREAST <- avg %>% 
  select(starts_with("BREAST"))

colnames(BREAST) <- BREAST[1,] 
BREAST <- BREAST[-1,]
BREAST <- as.data.frame(BREAST)
BREAST$malignant <- vi[,7]

UPPER_AERODIGESTIVE_TRACT <- avg %>% 
  select(starts_with("HEAD.AND.NECK"))

colnames(UPPER_AERODIGESTIVE_TRACT) <- UPPER_AERODIGESTIVE_TRACT[1,] 
UPPER_AERODIGESTIVE_TRACT <- UPPER_AERODIGESTIVE_TRACT[-1,]

STOMACH <- avg %>% 
  select(starts_with("STOMACH"))

colnames(STOMACH) <- STOMACH[1,] 
STOMACH <- STOMACH[-1,]

LARGE_INTESTINE <- avg %>% 
  select(starts_with("COLON.COLORECTAL"))

colnames(LARGE_INTESTINE) <- LARGE_INTESTINE[1,] 
LARGE_INTESTINE <- LARGE_INTESTINE[-1,]

PROSTATE <- avg %>% 
  select(starts_with("PROSTATE"))

colnames(PROSTATE) <- PROSTATE[1,] 
PROSTATE <- PROSTATE[-1,]

KIDNEY <- avg %>% 
  select(starts_with("KIDNEY"))

colnames(KIDNEY) <- KIDNEY[1,] 
KIDNEY <- KIDNEY[-1,]

LUNG <- avg %>% 
  select(starts_with("LUNG"))

colnames(LUNG) <- LUNG[1,] 
LUNG <- LUNG[-1,]

BONE <- avg %>% 
  select(starts_with("BONE"))

colnames(BONE) <- BONE[1,] 
BONE <- BONE[-1,]

PANCREAS <- avg %>% 
  select(starts_with("PANCREAS"))

colnames(PANCREAS) <- PANCREAS[1,] 
PANCREAS <- PANCREAS[-1,]
PANCREAS <- as.data.frame(PANCREAS)
PANCREAS$malignant <- vi[,1]


  OVARY <- avg %>% 
  select(starts_with("OVARY"))

colnames(OVARY) <- OVARY[1,] 
OVARY <- OVARY[-1,]

OESOPHAGUS <- avg %>% 
  select(starts_with("OESOPHAGUS"))

colnames(OESOPHAGUS) <- OESOPHAGUS[1,] 
OESOPHAGUS <- OESOPHAGUS[-1,]


PANCREAS <- mutate_all(PANCREAS, function(x) as.numeric(as.character(x)))



pancreas.malignant <- as.data.frame(cor(PANCREAS[, c(1,2,3,4,5,6)]))

rownames(pancreas.malignant) <- colnames(pancreas.malignant)

ggplot(pancreas.malignant, aes(x = rownames(pancreas.malignant), y = colnames(pancreas.malignant))) +
  geom_point() +
  labs(x = "Variable 1", y = "Variable 2", title = "Scatterplot with Correlation Line") +
  theme_bw()



ggplot(pancreas.malignant, aes(x = rownames(pancreas.malignant), y = colnames(pancreas.malignant), color = colnames(pancreas.malignant))) +
  geom_point() +  sm_statCorr()
  

ggplot(pancreas.malignant, aes(x = rownames(pancreas.malignant), y = pancreas.malignant$HUPT3)) +
  geom_point(aes(color = rownames(pancreas.malignant))) + sm_statCorr()


ggplot(pancreas.malignant, aes(x = rownames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$HUPT3, color = colnames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$HUPT4, color = colnames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$PANC0203, color = colnames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$ASPC1, color = colnames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$SU8686, color = colnames(pancreas.malignant))) +
  geom_point(aes(y = pancreas.malignant$malignant, color = colnames(pancreas.malignant))) +  
  labs(x = "Variable 1", y = "Variable 2", title = "Scatterplot with Correlation Line")


for (i in colnames(pancreas.malignant)) {
  print(i)
}


p1 <- ggscatter(PANCREAS, x = 'malignant', y = 'HUPT3', 
          add = "reg.line", cor.coef = TRUE, 
          cor.coef.label = paste0("correlation = ", round(pancreas.malignant, 2)))  
p2 <- ggscatter(PANCREAS, x = 'malignant', y = 'HUPT4', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(pancreas.malignant, 2))) 
p3 <- ggscatter(PANCREAS, x = 'malignant', y = 'PANC0203', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(pancreas.malignant, 2))) 
p4 <- ggscatter(PANCREAS, x = 'malignant', y = 'ASPC1', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(pancreas.malignant, 2))) 
p5 <- ggscatter(pancreas.malignant, x = 'malignant', y = 'SU8686', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(pancreas.malignant, 2))) 
p1 | p2 | p3 | p4 | p5

corrplot(cor(PANCREAS), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(PANCREAS), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(merged))

corrplot(cor(PANCREAS), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(PANCREAS), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(pancreas.malignant, c("everything", "pearson"), label = TRUE, insignificant = "blank")

BREAST <- mutate_all(BREAST, function(x) as.numeric(as.character(x)))


BREAST.malignant <- as.data.frame(cor(BREAST[, c(1,2,3,4,5,6)]))

p1 <- ggscatter(BREAST, x = 'malignant', y = 'BT474', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(BREAST.malignant, 2)))  
p2 <- ggscatter(BREAST, x = 'malignant', y = 'BT549', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(BREAST.malignant, 2))) 
p3 <- ggscatter(BREAST, x = 'malignant', y = 'EFM192A', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(BREAST.malignant, 2))) 
p4 <- ggscatter(BREAST, x = 'malignant', y = 'HCC1419', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(BREAST.malignant, 2))) 
p5 <- ggscatter(BREAST.malignant, x = 'malignant', y = 'CAMA1', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(BREAST.malignant, 2))) 
p1 | p2 | p3 | p4 | p5

corrplot(cor(BREAST), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(BREAST), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(merged))

corrplot(cor(BREAST), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(BREAST), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(BREAST.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

OVARY <- avg %>% 
  select(starts_with("OVARY"))

colnames(OVARY) <- OVARY[1,] 
OVARY <- OVARY[-1,]
OVARY$malignant <- vi[,3] 
OVARY <- OVARY[-39,] # Remover FGF21

OVARY <- mutate_all(OVARY, function(x) as.numeric(as.character(x)))

OVARY.malignant <- as.data.frame(cor(OVARY[, c(1,2,3,4,5,6)]))

OVARY.malignant

OVARY.malignant <- as.data.frame(cor(OVARY[, c(1,2,3,4,5,6)]))

p1 <- ggscatter(OVARY, x = 'malignant', y = 'OVSAHO', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(OVARY.malignant, 2)))  
p2 <- ggscatter(OVARY, x = 'malignant', y = 'JHOS2', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(OVARY.malignant, 2))) 
p3 <- ggscatter(OVARY, x = 'malignant', y = 'OAW28', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(OVARY.malignant, 2))) 
p4 <- ggscatter(OVARY, x = 'malignant', y = 'OVCAR4', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(OVARY.malignant, 2))) 
p5 <- ggscatter(OVARY, x = 'malignant', y = 'ONCODG1', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(OVARY.malignant, 2))) 
p1 | p2 | p3 

p4 | p5

corrplot(cor(OVARY), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(OVARY), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(OVARY.malignant))

corrplot(cor(OVARY), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(OVARY), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(OVARY.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

KIDNEY <- avg %>% 
  select(starts_with("KIDNEY"))

colnames(KIDNEY) <- KIDNEY[1,] 
KIDNEY <- KIDNEY[-1,]

KIDNEY$malignant <- vi[,4]

KIDNEY <- KIDNEY[-4,] # Remover ADIPOQ

KIDNEY <- mutate_all(KIDNEY, function(x) as.numeric(as.character(x)))

KIDNEY.malignant <- as.data.frame(cor(KIDNEY[, c(1,2,3,4,5,6)]))

KIDNEY.malignant

KIDNEY.malignant <- as.data.frame(cor(KIDNEY[, c(1,2,3,4,5,6)]))

p1 <- ggscatter(KIDNEY, x = 'malignant', y = 'RCC10RGB', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(KIDNEY.malignant, 2)))  
p2 <- ggscatter(KIDNEY, x = 'malignant', y = 'OSRC2', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(KIDNEY.malignant, 2))) 
p3 <- ggscatter(KIDNEY, x = 'malignant', y = 'VMRCRCZ', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(KIDNEY.malignant, 2))) 
p4 <- ggscatter(KIDNEY, x = 'malignant', y = 'KMRC3', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(KIDNEY.malignant, 2))) 
p5 <- ggscatter(KIDNEY, x = 'malignant', y = 'CAKI2', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(KIDNEY.malignant, 2))) 
p1 | p2 | p3 

p4 | p5

corrplot(cor(KIDNEY), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(KIDNEY), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(KIDNEY.malignant))

corrplot(cor(KIDNEY), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(KIDNEY), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(KIDNEY.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

SKIN <- avg %>% 
  select(starts_with("SKIN"))
colnames(SKIN) <- SKIN[1,] 
SKIN <- SKIN[-1,]
SKIN <- as.data.frame(SKIN)
SKIN$malignant <- vi[,8]
SKIN <- SKIN[c(-4, -32),] # Removendo ADIPOQ e CXCL8

SKIN <- mutate_all(SKIN, function(x) as.numeric(as.character(x)))

SKIN.malignant <- as.data.frame(cor(SKIN[, c(1,2,3,4,5,6)]))

SKIN.malignant

p1 <- ggscatter(SKIN, x = 'malignant', y = 'COLO741', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(SKIN.malignant, 2)))  
p2 <- ggscatter(SKIN, x = 'malignant', y = 'SKMEL3', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(SKIN.malignant, 2))) 
p3 <- ggscatter(SKIN, x = 'malignant', y = 'IGR1', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(SKIN.malignant, 2))) 
p4 <- ggscatter(SKIN, x = 'malignant', y = 'RVH421', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(SKIN.malignant, 2))) 
p5 <- ggscatter(SKIN, x = 'malignant', y = 'SKMEL5', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(SKIN.malignant, 2))) 
p1 | p2 | p3 
p4 | p5

corrplot(cor(SKIN), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(SKIN), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(SKIN.malignant))

corrplot(cor(SKIN), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(SKIN), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(SKIN.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

PROSTATE <- avg %>% 
  select(starts_with("PROSTATE"))

colnames(PROSTATE) <- PROSTATE[1,] 
PROSTATE <- PROSTATE[-1,]
PROSTATE$malignant <- vi[,6]

PROSTATE <- PROSTATE[c(-39, -84),] # Removendo FGF21 e NPPB

PROSTATE <- mutate_all(PROSTATE, function(x) as.numeric(as.character(x)))

PROSTATE.malignant <- as.data.frame(cor(PROSTATE[, c(1,2,3)]))

PROSTATE.malignant

corrplot(cor(PROSTATE), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(PROSTATE), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(PROSTATE.malignant))

corrplot(cor(PROSTATE), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(PROSTATE), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(PROSTATE.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

CENTRAL_NERVOUS_SYSTEM <- avg %>% 
  select(starts_with("BRAIN"))

colnames(CENTRAL_NERVOUS_SYSTEM) <- CENTRAL_NERVOUS_SYSTEM[1,] 
CENTRAL_NERVOUS_SYSTEM <- CENTRAL_NERVOUS_SYSTEM[-1,]

CENTRAL_NERVOUS_SYSTEM$malignant <- vi[,9] 

CENTRAL_NERVOUS_SYSTEM <- CENTRAL_NERVOUS_SYSTEM[-24,] # Removendo CSF2


CENTRAL_NERVOUS_SYSTEM <- mutate_all(CENTRAL_NERVOUS_SYSTEM, function(x) as.numeric(as.character(x)))

CENTRAL_NERVOUS_SYSTEM.malignant <- as.data.frame(cor(CENTRAL_NERVOUS_SYSTEM[, c(1,2,3,4,5,6)]))

CENTRAL_NERVOUS_SYSTEM.malignant

colnames(CENTRAL_NERVOUS_SYSTEM)[5]  <- 'aaa42MGBA'

p1 <- ggscatter(CENTRAL_NERVOUS_SYSTEM, x = 'malignant', y = 'KNS42', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(CENTRAL_NERVOUS_SYSTEM.malignant, 2)))  
p2 <- ggscatter(CENTRAL_NERVOUS_SYSTEM, x = 'malignant', y = 'TM31', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(CENTRAL_NERVOUS_SYSTEM.malignant, 2))) 
p3 <- ggscatter(CENTRAL_NERVOUS_SYSTEM, x = 'malignant', y = 'DKMG', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(CENTRAL_NERVOUS_SYSTEM.malignant, 2))) 
p4 <- ggscatter(CENTRAL_NERVOUS_SYSTEM, x = 'malignant', y = 'CCFSTTG1', 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(CENTRAL_NERVOUS_SYSTEM.malignant, 2))) 
p5 <- ggscatter(CENTRAL_NERVOUS_SYSTEM, x = 'malignant', y = "aaa42MGBA", 
                add = "reg.line", cor.coef = TRUE, 
                cor.coef.label = paste0("correlation = ", round(CENTRAL_NERVOUS_SYSTEM.malignant, 2))) 
p1 | p2 | p3 
p4 | p5

corrplot(cor(CENTRAL_NERVOUS_SYSTEM), method = 'shade', order = 'AOE', diag = FALSE)

corrplot(cor(CENTRAL_NERVOUS_SYSTEM), method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

corrplot(cor(CENTRAL_NERVOUS_SYSTEM.malignant))

corrplot(cor(CENTRAL_NERVOUS_SYSTEM), type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

corrplot(cor(CENTRAL_NERVOUS_SYSTEM), type="upper", order="hclust", diag = FALSE,
         col=brewer.pal(n=8, name="RdYlBu"))

ggcorr(CENTRAL_NERVOUS_SYSTEM.malignant, c("everything", "pearson"), label = TRUE,  insignificant = "blank")

ggcorr(cor(CENTRAL_NERVOUS_SYSTEM), c("everything", "pearson"), label = TRUE,  insignificant = "blank")      



ggplot(pancreas.malignant, aes(x = pancreas.malignant$malignant, y = c('HUPT3', 'HUPT4','PANC0203', 'ASPC1','SU8686'))) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
