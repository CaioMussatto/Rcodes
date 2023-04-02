remotes::install_github('JefersonSSouza/CaCao.seurat.plots')
library(seuratPercentages)


seurat_barplot(object = linhagens, # seu objeto seurat
               features_list = lista_CIFS[0:10], # sua lista de genes Ex. genes <- c('CCL8','IL6`,...)
               ident.use = unique(linhagens$sample, # ident a ser classificado nesse caso celltype
               show_percentage_legend =  T, # mostra porcentagem por ident TRUE ou FALSE
               percentage_legend_size= 2, # tamanho da fonte de porcentagem
               ncol = 8,  # n colunas do gráfico
               path_to_save = getwd(), # diretório a ser salvo o arquivo
               width = 18, # largura da imagem salva
               height = 25, #altura da imagem salva 
               plot_name = 'THYROID_barplot_features_1'))  # nome do arquivo a ser salvo


seurat_barplot(object = Lung, # seu objeto seurat
               features_list = lista_CIFS_lung$...1[41:87], # sua lista de genes Ex. genes <- c('CCL8','IL6`,...)
               ident.use = 'sample', # ident a ser classificado nesse caso celltype
               show_percentage_legend =  T, # mostra porcentagem por ident TRUE ou FALSE
               percentage_legend_size= 2, # tamanho da fonte de porcentagem
               ncol = 8,  # n colunas do gráfico
               path_to_save = getwd(), # diretório a ser salvo o arquivo
               width = 18, # largura da imagem salva
               height = 25, #altura da imagem salva 
               plot_name = 'Lung_barplot_features_2')



