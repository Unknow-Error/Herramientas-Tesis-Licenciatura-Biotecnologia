#Si no se tienen las librerías:

#install.packages("BiocManager")
#BiocManager::install(c("GenomicRanges", "rtracklayer", "ChIPseeker", "GenomicFeatures", "txdbmaker", "clusterProfiler"))
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("stringr")

#Invocar las librerías necesarias:
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(GenomicFeatures)
library(txdbmaker)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(stringr)

# Cargar FIMO en .tsv
#Factores de transcripción analizados: BbAP1, BbAreA, BbCdr1, BbCreA, BbCrz1, BbCtf1alpha, BbCtf1beta, BbFkh2, BbFlbC, BbFlbD, BbHcr1

factores_transcripcion <- c("BbAP1", "BbAreA", "BbCdr1", "BbCreA", "BbCrz1", "BbCtf1alpha", "BbCtf1beta",
                            "BbFkh2", "BbFlbC", "BbFlbD", "BbHcr1", "BbHox2", "BbHsf1", "BbKlf", "BbMb1",
                            "BbMcm1", "BbMsn2", "BbOps3", "BbOsrR1", "BbOsrR3", "BbOtf1", "BbPacC", "BbRei1",
                            "BbRep1", "BbSmr1", "BbSte12", "BbTenR", "BbYap1")

ruta_base_fimo <- "" #Indicar ruta de los archivos .tsv del FIMO

resultados_fimo <- list()

for (factor in factores_transcripcion){
  archivo <- paste0(ruta_base_fimo, "FIMO - ", factor, " - 1E-4.tsv")
  resultados_fimo[[factor]] <- read.table(archivo, header = TRUE, sep = "\t")
}

#Convertir data de FIMO en GenomicRanges

fimo_grs <- list()

for (factor in factores_transcripcion){
  fimo_grs[[factor]] <- GRanges(
    seqnames = resultados_fimo[[factor]]$sequence_name,
    ranges = IRanges(start = resultados_fimo[[factor]]$start, end = resultados_fimo[[factor]]$stop),
    strand = resultados_fimo[[factor]]$strand,
    score = resultados_fimo[[factor]]$score,
    p_value = resultados_fimo[[factor]][["p.value"]],
    q_value = resultados_fimo[[factor]][["q.value"]]
  )
}

#Cargar Anotaciones del genoma
gtf <- import("") #Indicar ruta y nombre del archivo .gtf que contiene las anotación del Genoma

#Anotar sitios de unión con ChipSeeker
txdb <- makeTxDbFromGFF(gtf, format = "gtf")

annotated_peaks <- list()

for (factor in factores_transcripcion){
  annotated_peaks[[factor]] <- annotatePeak(
    fimo_grs[[factor]],
    TxDb = txdb,
    tssRegion = c(-5000, 5000), #Cambiar a c(-3000, 3000) que es el estándar si se desea o a c(-1000, 1000) para analizar cajas estrictas del TSS
    level = "gene",
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
    sameStrand = FALSE,
    ignoreOverlap = FALSE,
    ignoreUpstream = FALSE,
    overlap = "TSS",
    verbose = TRUE,
    columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
}

# Convierte el resultado de annotatePeak a un DataFrame

anno_dfs <- list()

for (factor in factores_transcripcion){
  anno_dfs[[factor]] <- as.data.frame(annotated_peaks[[factor]])
}

#Chequear
head(anno_dfs[["BbAP1"]])
head(anno_dfs[["BbCtf1alpha"]])
head(anno_dfs[["BbTenR"]])

# Crear un gráfico de barras para comparar los factores
peak_count <- c()

for (factor in factores_transcripcion){
  peak_count <- c(peak_count, length(anno_dfs[[factor]]$annotation))
}

# Dataframe para el gráfico
df_peaks_FTs <- data.frame(Factor = factor(factores_transcripcion, levels = factores_transcripcion), Count = peak_count)

# Gráfico de barras mejorado
ggplot(df_peaks_FTs, aes(x = Factor, y = Count, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "plasma") +  # Paleta con más colores
  labs(x = "Factor de Transcripción", y = "Número de Picos",
       title = "Comparación de Cantidad de Picos FIMO entre Factores de Transcripción") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 25, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 25),
        axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 22, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank()) +
  coord_flip()  # Cambia la orientación del gráfico para mostrar mejor los


# Gráfico de caja para la puntuación FIMO de diferentes factores
#La puntuación FIMO indica qué tan fuerte es la coincidencia del motivo con la secuencia, por lo que comparar la distribución de las puntuaciones entre los FT y los controles puede ayudar a identificar diferencias en la afinidad del motivo entre los diferentes factores.

df_scores <- do.call(rbind, lapply(factores_transcripcion, function(factor) {
  data.frame(Factor = factor, Score = anno_dfs[[factor]][["score"]])
}))


ggplot(df_scores, aes(x = reorder(Factor, Score, FUN = median), y = Score)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "Factor de Transcripción (ordenado por mediana)", y = "Puntuación FIMO") +
  ggtitle("Comparación de Puntuaciones FIMO") +
  theme_classic() +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 36, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        legend.position = "none",  # Ocultar la leyenda ya que los nombres están en el eje X
        panel.grid.major.x = element_blank())
