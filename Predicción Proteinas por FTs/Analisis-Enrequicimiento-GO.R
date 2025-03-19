#BiocManager::install("topGO")

library(topGO)

#Análisis Enriquecimiento - GO


# Leer archivo GAF <- Contiene los GO annotations
gaf_file <- "" #Ruta del archivo .gaf con las anotaciones GO
gaf_data <- read.delim(gaf_file, header = FALSE, comment.char = "!", sep = "\t", stringsAsFactors = FALSE)

# Extraer columnas necesarias: ID del gen (DB_Object_ID) y términos GO (GO_ID)
colnames(gaf_data) <- c("DB", "GeneID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", "With_Or_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By")

# Verificar las primeras filas
head(gaf_data)

# Lista de genes anotados en los picos
# Suponiendo que tu archivo annotated-peaks tiene una columna `gene_id` con los IDs de NCBI
FT_gene_selection <- list()

for(factor in factores_transcripcion){
  gene_proteins_peaks <- read.table(paste(ruta_anotaciones_fimo, "anno.df_", factor, ".tsv", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  gene_list <- unique(gene_proteins_peaks[["GeneID"]])
  FT_gene_selection[[factor]] <- factor(as.integer(gaf_data$GeneID %in% gene_list))
}

for(factor in factores_transcripcion){
  names(FT_gene_selection[[factor]]) <- gaf_data$GeneID
}

# Crear un objeto topGOdata para análisis de enriquecimiento
GOData_FT_BP <- list()
GOData_FT_MF <- list()
GOData_FT_CC <- list()

for(factor in factores_transcripcion){
  GOData_FT_BP[[factor]] <- new("topGOdata",
    ontology = "BP",
    allGenes = FT_gene_selection[[factor]],
    nodeSize = 10, # nodeSize filtra términos GO con pocos genes
    annot = annFUN.gene2GO,
    gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
}

for(factor in factores_transcripcion){
  GOData_FT_MF[[factor]] <- new("topGOdata",
    ontology = "MF",
    allGenes = FT_gene_selection[[factor]],
    nodeSize = 10, # nodeSize filtra términos GO con pocos genes
    annot = annFUN.gene2GO,
    gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
}

for(factor in factores_transcripcion){
  GOData_FT_CC[[factor]] <- new("topGOdata",
    ontology = "CC",
    allGenes = FT_gene_selection[[factor]],
    nodeSize = 10, # nodeSize filtra términos GO con pocos genes
    annot = annFUN.gene2GO,
    gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
}

# Realizar el análisis usando el método 'classic' o 'elim' y el test estadístico Fisher
resultados_fisher_GO_FT_BP <- list()
resultados_fisher_GO_FT_MF <- list()
resultados_fisher_GO_FT_CC <- list()

for(factor in factores_transcripcion){
  resultados_fisher_GO_FT_BP[[factor]] <- runTest(GOData_FT_BP[[factor]], algorithm = "classic", statistic = "fisher")
}

for(factor in factores_transcripcion){
  resultados_fisher_GO_FT_MF[[factor]] <- runTest(GOData_FT_MF[[factor]], algorithm = "classic", statistic = "fisher")
}

for(factor in factores_transcripcion){
  resultados_fisher_GO_FT_CC[[factor]] <- runTest(GOData_FT_CC[[factor]], algorithm = "classic", statistic = "fisher")
}

top_GO_results_FT_BP <- list()
top_GO_results_FT_MF <- list()
top_GO_results_FT_CC <- list()

for(factor in factores_transcripcion){
  top_GO_results_FT_BP[[factor]] <- GenTable(GOData_FT_BP[[factor]],
    classicFisher = resultados_fisher_GO_FT_BP[[factor]],
    orderBy = "classicFisher",
    topNodes = 6)
}

for(factor in factores_transcripcion){
  top_GO_results_FT_MF[[factor]] <- GenTable(GOData_FT_MF[[factor]],
    classicFisher = resultados_fisher_GO_FT_MF[[factor]],
    orderBy = "classicFisher",
    topNodes = 6)
}

for(factor in factores_transcripcion){
  top_GO_results_FT_CC[[factor]] <- GenTable(GOData_FT_CC[[factor]],
    classicFisher = resultados_fisher_GO_FT_CC[[factor]],
    orderBy = "classicFisher",
    topNodes = 6)
}

# Crear un gráfico de barras con los términos GO enriquecidos
ruta_graficos_GO <- ""

# Combinar los resultados de BP en un solo dataframe
factores_seleccionados <- c("BbAP1", "BbCreA", "BbCrz1", "BbFlbC", "BbHcr1", "BbKlf", "BbOps3", "BbOsrR1", "BbOsrR3", "BbOtf1",  "BbRei1", "BbRep1")
factores_seleccionados_2 <- c("BbCtf1alpha","BbHox2", "BbSmr1", "BbTenR", "BbAreA")
factores_seleccionados_3 <- c("BbCdr1","BbCtf1beta", "BbFkh2", "BbFlbD","BbYap1")
factores_seleccionados_4 <- c ("BbHsf1", "BbMb1", "BbMsn2", "BbPacC","BbSte12" )

df_BP_Select_1 <- do.call(rbind, lapply(factores_seleccionados, function(factor) {
  df <- top_GO_results_FT_BP[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_BP_Select_2 <- do.call(rbind, lapply(factores_seleccionados_2, function(factor) {
  df <- top_GO_results_FT_BP[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_BP_Select_3 <- do.call(rbind, lapply(factores_seleccionados_3, function(factor) {
  df <- top_GO_results_FT_BP[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_BP_Select_4 <- do.call(rbind, lapply(factores_seleccionados_4, function(factor) {
  df <- top_GO_results_FT_BP[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

# Convertir p-valor a -log10(p-value) para mejor visualización
df_BP_Select_1$logP <- -log10(as.numeric(df_BP_Select_1$classicFisher))
df_BP_Select_2$logP <- -log10(as.numeric(df_BP_Select_2$classicFisher))
df_BP_Select_3$logP <- -log10(as.numeric(df_BP_Select_3$classicFisher))
df_BP_Select_4$logP <- -log10(as.numeric(df_BP_Select_4$classicFisher))

# Graficar Dot Plot
ggplot(df_BP_Select_1, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Procesos Biológicos)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 40, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_BP_Select_1_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_BP_Select_2, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Procesos Biológicos)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_BP_Select_2_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_BP_Select_3, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Procesos Biológicos)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_BP_Select_3_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_BP_Select_4, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Procesos Biológicos)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_BP_Select_4_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)


df_MF_Select_1 <- do.call(rbind, lapply(factores_seleccionados, function(factor) {
  df <- top_GO_results_FT_MF[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_MF_Select_2 <- do.call(rbind, lapply(factores_seleccionados_2, function(factor) {
  df <- top_GO_results_FT_MF[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_MF_Select_3 <- do.call(rbind, lapply(factores_seleccionados_3, function(factor) {
  df <- top_GO_results_FT_MF[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_MF_Select_4 <- do.call(rbind, lapply(factores_seleccionados_4, function(factor) {
  df <- top_GO_results_FT_MF[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

# Convertir p-valor a -log10(p-value) para mejor visualización
df_MF_Select_1$logP <- -log10(as.numeric(df_MF_Select_1$classicFisher))
df_MF_Select_2$logP <- -log10(as.numeric(df_MF_Select_2$classicFisher))
df_MF_Select_3$logP <- -log10(as.numeric(df_MF_Select_3$classicFisher))
df_MF_Select_4$logP <- -log10(as.numeric(df_MF_Select_4$classicFisher))

# Graficar Dot Plot
ggplot(df_MF_Select_1, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Función Molecular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 40, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_MF_Select_1_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_MF_Select_2, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Función Molecular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_MF_Select_2_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_MF_Select_3, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Función Molecular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_MF_Select_3_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_MF_Select_4, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Función Molecular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_MF_Select_4_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

df_CC_Select_1 <- do.call(rbind, lapply(factores_seleccionados, function(factor) {
  df <- top_GO_results_FT_CC[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_CC_Select_2 <- do.call(rbind, lapply(factores_seleccionados_2, function(factor) {
  df <- top_GO_results_FT_CC[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_CC_Select_3 <- do.call(rbind, lapply(factores_seleccionados_3, function(factor) {
  df <- top_GO_results_FT_CC[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

df_CC_Select_4 <- do.call(rbind, lapply(factores_seleccionados_3, function(factor) {
  df <- top_GO_results_FT_CC[[factor]]
  df$Factor <- factor  # Agregar la columna de factor de transcripción
  return(df)
}))

# Convertir p-valor a -log10(p-value) para mejor visualización
df_CC_Select_1$logP <- -log10(as.numeric(df_CC_Select_1$classicFisher))
df_CC_Select_2$logP <- -log10(as.numeric(df_CC_Select_2$classicFisher))
df_CC_Select_3$logP <- -log10(as.numeric(df_CC_Select_3$classicFisher))
df_CC_Select_4$logP <- -log10(as.numeric(df_CC_Select_4$classicFisher))

# Graficar Dot Plot
ggplot(df_CC_Select_1, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Localización celular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 40, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_CC_Select_1_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_CC_Select_2, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Localización celular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_CC_Select_2_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_CC_Select_3, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Localización celular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_CC_Select_3_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)

ggplot(df_CC_Select_4, aes(x = Factor, y = Term, size = Significant, color = logP)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "Análisis Enriquecimiento GO (Localización celular)",
       x = "Factor de Transcripción",
       y = "Término GO",
       size = "Genes",
       color = "-log10(p-valor)") +
  theme_bw() +
  theme(
      plot.title = element_text(hjust = 0.5, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D", angle = 90),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 27, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 20, 10, 20),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))
ggsave(paste(ruta_graficos_GO, "GO_CC_Select_4_DotPlot_Top6.png", sep = ""), width = 24, height = 16, dpi = 300)
