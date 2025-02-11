library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(GenomicFeatures)
library(txdbmaker)
library(clusterProfiler)
library(ReactomePA)
library(topGO)
library(biomaRt)
library(ggplot2)
library(Rgraphviz)
library(dplyr)
library(stringr)
library(VennDiagram)
library(RColorBrewer)
library(readr)
library(Biostrings)  # Para manejar archivos FASTA

#Convertir de PWM de CIS-BP a .meme
convert_pwm_to_meme <- function(pwm_file, output_meme, motif_name="TF_NAME") {
  # Cargar la matriz PWM desde el archivo (asumiendo formato CIS-BP)
  pwm_data <- read.delim(pwm_file, header=FALSE, row.names=1, sep="\t")

  # Obtener el número de posiciones (w)
  w <- ncol(pwm_data)

  # Crear el archivo .meme
  cat("MEME version 4\n\n", file=output_meme)
  cat("ALPHABET= ACGT\n\n", file=output_meme, append=TRUE)
  cat("strands: + -\n\n", file=output_meme, append=TRUE)
  cat(paste("MOTIF", motif_name, "\n\n"), file=output_meme, append=TRUE)
  cat(paste("letter-probability matrix: alength= 4 w=", w, "\n"), file=output_meme, append=TRUE)

  # Guardar la matriz PWM en formato MEME
  write.table(t(pwm_data), file=output_meme, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

# Uso de la función
convert_pwm_to_meme("BbHCR1-Neurospora crassa -NCU06407-M02614_2.00.tsv", "BbHCR1-Neurospora crassa -NCU06407-M02614_2.00.meme", motif_name="BbHCR1")


# Cargar FIMO en .tsv
#Factores de transcripción analizados: BbMSN2, BbNDT80, BbCreA, BbFOXK2, BbOtf1, BbHCR1, BbPacC
#Factores de transcripción Housekeeping (Control) : BbTBP.

fimo_results_BbMSN2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-BbMsn2-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")
fimo_results_BbNDT80 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO - BbNdt80-PhoG-Neurospora crassa-M02387_2.00.tsv", header = TRUE, sep = "\t")
fimo_results_BbFOXK2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-BbFOXK2-N-crassa-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")
fimo_results_BbPacC <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-PacC-N-crassa-Genomic ASM28067v1-1E-5.tsv", header = TRUE, sep = "\t")
fimo_results_BbOtf1 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-BbOtf1 - N.crassa-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")
fimo_results_BbHCR1 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-BbHCR1-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")

fimo_results_BbCreA <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-BbCreA-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")
fimo_results_BbTBP <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/FIMO-TATABP-Genomic ASM28067v1-1E-4.tsv", header = TRUE, sep = "\t")

#Convertir data de FIMO en GenomicRanges

fimo_gr_BbMSN2 <- GRanges(
  seqnames = fimo_results_BbMSN2$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbMSN2$start, end = fimo_results_BbMSN2$stop),
  strand = fimo_results_BbMSN2$strand,
  score = fimo_results_BbMSN2$score,
  p_value = fimo_results_BbMSN2[["p.value"]],
  q_value = fimo_results_BbMSN2[["q.value"]]
)

fimo_gr_BbNDT80 <- GRanges(
  seqnames = fimo_results_BbNDT80$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbNDT80$start, end = fimo_results_BbNDT80$stop),
  strand = fimo_results_BbNDT80$strand,
  score = fimo_results_BbNDT80$score,
  p_value = fimo_results_BbNDT80[["p.value"]],
  q_value = fimo_results_BbNDT80[["q.value"]]
)

fimo_gr_BbCreA <- GRanges(
  seqnames = fimo_results_BbCreA$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbCreA$start, end = fimo_results_BbCreA$stop),
  strand = fimo_results_BbCreA$strand,
  score = fimo_results_BbCreA$score,
  p_value = fimo_results_BbCreA[["p.value"]],
  q_value = fimo_results_BbCreA[["q.value"]]
)

fimo_gr_BbFOXK2 <- GRanges(
  seqnames = fimo_results_BbFOXK2$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbFOXK2$start, end = fimo_results_BbFOXK2$stop),
  strand = fimo_results_BbFOXK2$strand,
  score = fimo_results_BbFOXK2$score,
  p_value = fimo_results_BbFOXK2[["p.value"]],
  q_value = fimo_results_BbFOXK2[["q.value"]]
)

fimo_gr_BbPacC <- GRanges(
  seqnames = fimo_results_BbPacC$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbPacC$start, end = fimo_results_BbPacC$stop),
  strand = fimo_results_BbPacC$strand,
  score = fimo_results_BbPacC$score,
  p_value = fimo_results_BbPacC[["p.value"]],
  q_value = fimo_results_BbPacC[["q.value"]]
)

fimo_gr_BbOtf1 <- GRanges(
  seqnames = fimo_results_BbOtf1$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbOtf1$start, end = fimo_results_BbOtf1$stop),
  strand = fimo_results_BbOtf1$strand,
  score = fimo_results_BbOtf1$score,
  p_value = fimo_results_BbOtf1[["p.value"]],
  q_value = fimo_results_BbOtf1[["q.value"]]
)

fimo_gr_BbHCR1 <- GRanges(
  seqnames = fimo_results_BbHCR1$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbHCR1$start, end = fimo_results_BbHCR1$stop),
  strand = fimo_results_BbHCR1$strand,
  score = fimo_results_BbHCR1$score,
  p_value = fimo_results_BbHCR1[["p.value"]],
  q_value = fimo_results_BbHCR1[["q.value"]]
)

fimo_gr_BbTBP <- GRanges(
  seqnames = fimo_results_BbTBP$sequence_name,  # or appropriate column for chromosome
  ranges = IRanges(start = fimo_results_BbTBP$start, end = fimo_results_BbTBP$stop),
  strand = fimo_results_BbTBP$strand,
  score = fimo_results_BbTBP$score,
  p_value = fimo_results_BbTBP[["p.value"]],
  q_value = fimo_results_BbTBP[["q.value"]]
)

#Cargar Anotaciones del genoma
gtf <- import("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gtf")

#Anotar sitios de unión con ChipSeeker
txdb <- makeTxDbFromGFF("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gtf", format = "gtf")
annotated_peaks_BbOtf1 <- annotatePeak(fimo_gr_BbOtf1,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbMSN2 <- annotatePeak(fimo_gr_BbMSN2,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbNDT80 <- annotatePeak(fimo_gr_BbNDT80,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbHCR1 <- annotatePeak(fimo_gr_BbHCR1,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbCreA <- annotatePeak(fimo_gr_BbCreA,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbFOXK2 <- annotatePeak(fimo_gr_BbFOXK2,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbPacC <- annotatePeak(fimo_gr_BbPacC,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))
annotated_peaks_BbTBP <- annotatePeak(fimo_gr_BbTBP,
                                      TxDb = txdb,
                                      tssRegion = c(-1000, 1000),
                                      level = "gene",
                                      assignGenomicAnnotation = TRUE,
                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                                      sameStrand = FALSE,
                                      ignoreOverlap = FALSE,
                                      ignoreUpstream = FALSE,
                                      overlap = "TSS",
                                      verbose = TRUE,
                                      columns = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"))

# Convierte el resultado de annotatePeak a un DataFrame
anno.df_BbMSN2 <- as.data.frame(annotated_peaks_BbMSN2)
anno.df_BbNDT80 <- as.data.frame(annotated_peaks_BbNDT80)
anno.df_BbHCR1 <- as.data.frame(annotated_peaks_BbHCR1)
anno.df_BbCreA <- as.data.frame(annotated_peaks_BbCreA)
anno.df_BbFOXK2 <- as.data.frame(annotated_peaks_BbFOXK2)
anno.df_BbTBP <- as.data.frame(annotated_peaks_BbTBP)
anno.df_BbPacC <- as.data.frame(annotated_peaks_BbPacC)
anno.df_BbOtf1 <- as.data.frame(annotated_peaks_BbOtf1)

#Inspeccionar y visualizar
head(anno.df_BbMSN2)
head(anno.df_BbNDT80)
head(anno.df_BbHCR1)
head(anno.df_BbCreA)
head(anno.df_BbFOXK2)
head(anno.df_BbTBP)
head(anno.df_BbPacC)
head(anno.df_BbOtf1)

# Crear un gráfico de barras para comparar los factores
factor_list <- c("BbMSN2", "BbNDT80", "BbHCR1", "BbFOXK2", "BbPacC", "BbOtf1", "BbCreA", "BbTBP")
peak_count <- c(length(anno.df_BbMSN2$annotation), length(anno.df_BbNDT80$annotation), length(anno.df_BbHCR1$annotation),
                length(anno.df_BbFOXK2$annotation), length(anno.df_BbPacC$annotation), length(anno.df_BbOtf1$annotation),
                length(anno.df_BbCreA$annotation), length(anno.df_BbTBP$annotation))

# Dataframe para el gráfico
df_peaks_FTs <- data.frame(Factor = factor(factor_list, levels = factor_list), Count = peak_count)

# Gráfico de barras
ggplot(df_peaks_FTs, aes(x = Factor, y = Count, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Factor de Transcripción", y = "Número de Picos",
       title = "Comparación de Cantidad de Picos entre Factores de Transcripción") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 36, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        legend.position = "none",  # Ocultar la leyenda ya que los nombres están en el eje X
        panel.grid.major.x = element_blank())

#Mismo gráfico pero filtrando solo los "picos de los promotores"
promoter_peaks_BbMSN2 <- subset(anno.df_BbMSN2, annotation == "Promoter")
promoter_peaks_BbNDT80 <- subset(anno.df_BbNDT80, annotation == "Promoter")
promoter_peaks_BbHCR1 <- subset(anno.df_BbHCR1, annotation == "Promoter")
promoter_peaks_BbFOXK2 <- subset(anno.df_BbFOXK2, annotation == "Promoter")
promoter_peaks_BbPacC <- subset(anno.df_BbPacC, annotation == "Promoter")
promoter_peaks_BbOtf1 <- subset(anno.df_BbOtf1, annotation == "Promoter")
promoter_peaks_BbCreA <- subset(anno.df_BbCreA, annotation == "Promoter")
promoter_peaks_BbTBP <- subset(anno.df_BbTBP, annotation == "Promoter")

promoter_peak_count <- c(length(promoter_peaks_BbMSN2$annotation), length(promoter_peaks_BbNDT80$annotation),
                         length(promoter_peaks_BbHCR1$annotation), length(promoter_peaks_BbFOXK2$annotation), length(promoter_peaks_BbPacC$annotation),
                         length(promoter_peaks_BbOtf1$annotation), length(promoter_peaks_BbCreA$annotation), length(promoter_peaks_BbTBP$annotation))

df_promotor_peaks_FTs <- data.frame(Factor = factor(factor_list, levels = factor_list), Count = promoter_peak_count)

ggplot(df_promotor_peaks_FTs, aes(x = Factor, y = Count, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Factor de Transcripción", y = "Número de Promotores",
       title = "Comparación de Cantidad de Promotores activados según Factor de Transcripción") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 36, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        legend.position = "none",  # Ocultar la leyenda ya que los nombres están en el eje X
        panel.grid.major.x = element_blank())

# Gráfico de caja para la puntuación FIMO de diferentes factores
#La puntuación FIMO indica qué tan fuerte es la coincidencia del motivo con la secuencia, por lo que comparar la distribución de las puntuaciones entre los FT y los controles puede ayudar a identificar diferencias en la afinidad del motivo entre los diferentes factores.
df_scores <- rbind(
  data.frame(Factor = "BbMSN2", Score = anno.df_BbMSN2[["score"]]),
  data.frame(Factor = "BbNDT80", Score = anno.df_BbNDT80[["score"]]),
  data.frame(Factor = "BbHCR1", Score = anno.df_BbHCR1[["score"]]),
  data.frame(Factor = "BbFOXK2", Score = anno.df_BbFOXK2[["score"]]),
  data.frame(Factor = "BbPacC", Score = anno.df_BbPacC[["score"]]),
  data.frame(Factor = "BbOtf1", Score = anno.df_BbOtf1[["score"]]),
  data.frame(Factor = "BbCreA", Score = anno.df_BbCreA[["score"]]),
  data.frame(Factor = "BbTBP", Score = anno.df_BbTBP[["score"]])
)

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


# Diagrama de Venn para ver la co-ocurrencia de picos entre dos factores
# Un diagrama de Venn te puede ayudar a visualizar qué factores de transcripción comparten picos en las mismas ubicaciones.
venn_list_1 <- list(
  BbMSN2 = unique(anno.df_BbMSN2$geneId),
  BbNDT80 = unique(anno.df_BbNDT80$geneId),
  BbHCR1 = unique(anno.df_BbHCR1$geneId),
  BbFOXK2 = unique(anno.df_BbFOXK2$geneId),
  BbOtf1 = unique(anno.df_BbOtf1$geneId)
)

myCol <- brewer.pal(5, "Pastel2")
venn.diagram(venn_list_1, main = "Co-ocurrencia de Picos entre BbMSN2, BbNDT80, BbHCR1, BbFOXK2  y BbOtf1",
            output = TRUE,
            filename = "Diagrama-Venn-FTs-1.png",
            imagetype = "png",
            lwd = 2,
            col = myCol,
            fill = myCol,
            fontface = 'bold',
            fontfamily = 'sans',
            cat.cex = '0.6',
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.fontfamily = "sans")

venn_list_2 <- list(
  BbCreA = unique(anno.df_BbCreA$geneId),
  BbTBP = unique(anno.df_BbTBP$geneId)
)

myCol_2 <- brewer.pal(3, "Pastel2")[1:2]
venn.diagram(venn_list_2, main = "Co-ocurrencia de Picos entre los controles BbCreA y BbTBP",
            output = TRUE,
            filename = "Diagrama-Venn-FTs-2.png",
            imagetype = "png",
            lwd = 2,
            col = myCol_2,
            fill = myCol_2,
            fontface = 'bold',
            fontfamily = 'sans',
            cat.cex = '0.6',
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.fontfamily = "sans")

venn_list_3 <- list(
  BbNDT80 = unique(anno.df_BbNDT80$geneId),
  BbHCR1 = unique(anno.df_BbHCR1$geneId),
  BbFOXK2 = unique(anno.df_BbFOXK2$geneId),
  BbOtf1 = unique(anno.df_BbOtf1$geneId),
  BbPacC = unique(anno.df_BbPacC$geneId)
)

venn.diagram(venn_list_3, main = "Co-ocurrencia de Picos entre BbNDT80, BbHCR1, BbFOXK2, BbPacC  y BbOtf1",
            output = TRUE,
            filename = "Diagrama-Venn-FTs-3.png",
            imagetype = "png",
            lwd = 2,
            col = myCol,
            fill = myCol,
            fontface = 'bold',
            fontfamily = 'sans',
            cat.cex = '0.6',
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.fontfamily = "sans")

library(ggVennDiagram)

venn_list_4 <- list(
  BbMSN2 = unique(anno.df_BbMSN2$geneId),
  BbNDT80 = unique(anno.df_BbNDT80$geneId),
  BbHCR1 = unique(anno.df_BbHCR1$geneId),
  BbFOXK2 = unique(anno.df_BbFOXK2$geneId),
  BbOtf1 = unique(anno.df_BbOtf1$geneId),
  BbPacC = unique(anno.df_BbPacC$geneId)
)

# Crear el diagrama de Venn
ggVennDiagram(venn_list_4) +
  ggtitle("Co-ocurrencia de Picos entre los FTs seleccionados") +
  theme_minimal()



# Gráfico de densidad de la distancia a los TSS
# Compara la localización de los picos respecto a los TSS para los diferentes factores de transcripción y su relación con los controles (TATABP y RAP1). Esto te ayudará a observar si hay una tendencia de los FT a ubicarse en regiones más cercanas a los TSS.
df_TSS <- rbind(
  data.frame(Factor = "BbMSN2", distanceToTSS = anno.df_BbMSN2[["distanceToTSS"]]),
  data.frame(Factor = "BbNDT80", distanceToTSS = anno.df_BbNDT80[["distanceToTSS"]]),
  data.frame(Factor = "BbHCR1", distanceToTSS = anno.df_BbHCR1[["distanceToTSS"]]),
  data.frame(Factor = "BbFOXK2", distanceToTSS = anno.df_BbFOXK2[["distanceToTSS"]]),
  data.frame(Factor = "BbPacC", distanceToTSS = anno.df_BbPacC[["distanceToTSS"]]),
  data.frame(Factor = "BbOtf1", distanceToTSS = anno.df_BbOtf1[["distanceToTSS"]]),
  data.frame(Factor = "BbCreA", distanceToTSS = anno.df_BbCreA[["distanceToTSS"]]),
  data.frame(Factor = "BbTBP", distanceToTSS = anno.df_BbTBP[["distanceToTSS"]])
)

ggplot(df_TSS, aes(x = distanceToTSS)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  facet_wrap(~ Factor, ncol = 4, scales = "free") +
  labs(x = "Distancia al TSS (pb)", y = "Densidad",
       title = "Distribución de la distancia al TSS por factor de transcripción") +
  theme(plot.title = element_text(hjust = 0.5, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        legend.position = "none",  # Ocultar la leyenda ya que los nombres están en el eje X
        panel.grid.major.x = element_blank())


ggplot(df_TSS, aes(x = Factor, y = distanceToTSS, fill = Factor)) +
  geom_violin() +
  labs(x = "Factor de Transcripción", y = "Distancia al TSS (pb)", title = "Distribución de la distancia al TSS por factor de transcripción") +
  theme(plot.title = element_text(hjust = 0.5, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        legend.position = "none",  # Ocultar la leyenda ya que los nombres están en el eje X
        panel.grid.major.x = element_blank())

#Filtrar y obtener las proteínas que transcriben los Factores de Transcrición en Beauveria bassiana

#Transformar GeneIDs en locus_Tag en .GAF file

# Definir la función
extract_geneid_locus_tag <- function(gbff_file, output_file) {
  # Leer el archivo GBFF como texto
  gbff_lines <- readLines(gbff_file)

  # Inicializar listas para almacenar locus_tags y GeneIDs
  locus_tags <- c()
  gene_ids <- c()

  # Recorrer las líneas en el archivo GBFF
  for (line in gbff_lines) {
    if (grepl("/locus_tag=", line)) {
      locus_tag <- str_extract(line, "(?<=/locus_tag=\")[^\"]+")
      locus_tags <- c(locus_tags, locus_tag)
    }
    if (grepl("/db_xref=\"GeneID:", line)) {
      gene_id <- str_extract(line, "(?<=GeneID:)[0-9]+")
      gene_ids <- c(gene_ids, gene_id)
    }
  }

  # Combinar en un data frame
  gene_mapping <- data.frame(GeneID = gene_ids, LocusTag = locus_tags, stringsAsFactors = FALSE)

  # Eliminar duplicados
  gene_mapping_unique <- unique(gene_mapping)

  # Guardar el resultado en un archivo .tsv
  write.table(gene_mapping_unique, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Obtener correlacion GeneIds <-> locus_tag
extract_geneid_locus_tag("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic_genID_LocusTag.tsv")

# Obtener un archivo que correlacione el locus_tag con GeneId, protein_id, transcript_id y product.
extract_gene_info_from_gtf <- function(gtf_file, output_file) {
  # Leer el archivo GTF
  gtf_lines <- readLines(gtf_file)

  # Inicializar una lista vacía para almacenar los resultados
  results <- list()

  # Recorrer cada línea del archivo
  for (line in gtf_lines) {
    if (grepl("^#", line) | line == "") next  # Ignorar comentarios y líneas vacías

    # Extraer la columna de atributos (última columna)
    attributes <- str_split(line, "\t")[[1]][9]

    # Usar expresiones regulares más específicas para evitar "orig_protein_id"
    locus_tag <- str_extract(attributes, '\\blocus_tag "([^"]+)"')
    gene_id <- str_extract(attributes, '\\bdb_xref "GeneID:([0-9]+)"')
    transcript_id <- str_extract(attributes, '\\btranscript_id "([^"]+)"')
    protein_id <- str_extract(attributes, '\\bprotein_id "([^"]+)"(?![a-zA-Z])')  # Asegura que no esté seguido de letras
    product <- str_extract(attributes, '\\bproduct "([^"]+)"')

    # Limpiar comillas y formato de las capturas
    locus_tag <- str_remove_all(locus_tag, 'locus_tag "|\"')
    gene_id <- str_remove_all(gene_id, 'db_xref "GeneID:|\"')
    transcript_id <- str_remove_all(transcript_id, 'transcript_id "|\"')
    protein_id <- str_remove_all(protein_id, 'protein_id "|\"')
    product <- str_remove_all(product, 'product "|\"')

    # Almacenar solo si transcript_id, protein_id y product no son NA
    if (!is.na(locus_tag) & !is.na(gene_id) &
        !is.na(transcript_id) & !is.na(protein_id) & !is.na(product)) {
      results <- append(results, list(c(locus_tag, gene_id, transcript_id, protein_id, product)))
    }
  }

  # Convertir la lista en un data frame
  gene_info <- do.call(rbind, results)
  gene_info_unique <- unique(gene_info)
  colnames(gene_info_unique) <- c("locus_tag", "GeneID", "transcript_id", "protein_id", "product")

  # Guardar el resultado en un archivo .tsv
  write.table(gene_info_unique, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Archivo generado con éxito:", output_file, "\n")
}

extract_gene_info_from_gtf("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gtf", "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic_gene_info.tsv")

gene_mapping <- read.table("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic_gene_info.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_mapping_2 <- read.table("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic_genID_LocusTag.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

head(gene_mapping)
head(gene_mapping_2)

anno.df_BbMSN2_2 <- merge(anno.df_BbMSN2, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbMSN2_2 <- anno.df_BbMSN2_2 %>% rename(locus_tag = geneId)
anno.df_BbNDT80_2 <- merge(anno.df_BbNDT80, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbNDT80_2 <- anno.df_BbNDT80_2 %>% rename(locus_tag = geneId)
anno.df_BbHCR1_2 <- merge(anno.df_BbHCR1, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbHCR1_2 <- anno.df_BbHCR1_2 %>% rename(locus_tag = geneId)
anno.df_BbFOXK2_2 <- merge(anno.df_BbFOXK2, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbFOXK2_2 <- anno.df_BbFOXK2_2 %>% rename(locus_tag = geneId)
anno.df_BbPacC_2 <- merge(anno.df_BbPacC, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbPacC_2 <- anno.df_BbPacC_2 %>% rename(locus_tag = geneId)
anno.df_BbOtf1_2 <- merge(anno.df_BbOtf1, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbOtf1_2 <- anno.df_BbOtf1_2 %>% rename(locus_tag = geneId)
anno.df_BbCreA_2 <- merge(anno.df_BbCreA, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbCreA_2 <- anno.df_BbCreA_2 %>% rename(locus_tag = geneId)
anno.df_BbTBP_2  <- merge(anno.df_BbTBP, gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
anno.df_BbTBP_2 <- anno.df_BbTBP_2 %>% rename(locus_tag = geneId)

write.table(anno.df_BbMSN2_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbMSN2.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbMSN2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbMSN2_origin.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbNDT80_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbNDT80.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbNDT80, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbNDT80_origin.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbHCR1_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbHCR1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbHCR1, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbHCR1_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbFOXK2_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbFOXK2.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbFOXK2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbFOXK2_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbPacC_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbPacC.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbPacC, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbPacC_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbOtf1_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbOtf1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbOtf1, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbOtf1_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbCreA_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbCreA.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbCreA, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbCreA_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbTBP_2, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbTBP.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(anno.df_BbTBP, "/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbTBP_original.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Comprobar lo guardado
anno.df_BbMSN2_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbMSN2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbNDT80_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbNDT80.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbPacC_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbPacC.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbFOXK2_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbFOXK2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbOtf1_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbOtf1.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbHCR1_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbHCR1.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbTBP_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbTBP.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anno.df_BbCreA_2 <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbCreA.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Extraer los genes proteicos
genesIDs_BbMSN2 <- anno.df_BbMSN2_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbMSN2 <- na.omit(genesIDs_BbMSN2)
genesIDs_BbMSN2 <- unique(genesIDs_BbMSN2)
genesIDs_BbNDT80 <- anno.df_BbNDT80_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbNDT80 <- na.omit(genesIDs_BbNDT80)
genesIDs_BbNDT80 <- unique(genesIDs_BbNDT80)
genesIDs_BbPacC <- anno.df_BbPacC_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbPacC <- na.omit(genesIDs_BbPacC)
genesIDs_BbPacC <- unique(genesIDs_BbPacC)
genesIDs_BbFOXK2 <- anno.df_BbFOXK2_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbFOXK2 <- na.omit(genesIDs_BbFOXK2)
genesIDs_BbFOXK2 <- unique(genesIDs_BbFOXK2)
genesIDs_BbOtf1 <- anno.df_BbOtf1_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbOtf1 <- na.omit(genesIDs_BbOtf1)
genesIDs_BbOtf1 <- unique(genesIDs_BbOtf1)
genesIDs_BbHCR1 <- anno.df_BbHCR1_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbHCR1 <- na.omit(genesIDs_BbHCR1)
genesIDs_BbHCR1 <- unique(genesIDs_BbHCR1)
genesIDs_BbTBP <- anno.df_BbTBP_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbTBP <- na.omit(genesIDs_BbTBP)
genesIDs_BbTBP <- unique(genesIDs_BbTBP)
genesIDs_BbCreA <- anno.df_BbCreA_2 %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
genesIDs_BbCreA <- na.omit(genesIDs_BbCreA)
genesIDs_BbCreA <- unique(genesIDs_BbCreA)

#Generar un .tsv de cada Factor de transcripción que contenga las proteínas las cuales induce su transcripción.
write.table(genesIDs_BbMSN2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbMSN2_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbNDT80, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbNDT80_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbPacC, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbPacC_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbFOXK2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbFOXK2_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbOtf1, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbOtf1_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbHCR1, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbHCR1_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbTBP, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbTBP_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(genesIDs_BbCreA, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbCreA_proteinas.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Otra función para generar un .csv que tenga los annotation peaks con los valores de annotacion protein_id, transcript_id, locus_tag, etc.
extract_annotations <- function(data_frame, gbff_path, output_csv) {
  
  # Read the .gbff file
  gbff_content <- readLines(gbff_path)
  
  # Parse annotations from the .gbff file
  annotation_list <- list()
  gene_id <- locus_tag <- transcript_id <- protein_id <- NA
  
  for (line in gbff_content) {
    
    # Check for GeneID
    if (grepl("/db_xref=\"GeneID:", line)) {
      gene_id <- sub(".*/db_xref=\"GeneID:(\\d+)\".*", "\\1", line)
    }
    
    # Check for locus_tag
    if (grepl("/locus_tag=", line)) {
      locus_tag <- sub(".*/locus_tag=\"([^\"]+)\".*", "\\1", line)
    }
    
    # Check for transcript_id
    if (grepl("/transcript_id=", line)) {
      transcript_id <- sub(".*/transcript_id=\"([^\"]+)\".*", "\\1", line)
    }
    
    # Check for protein_id if it exists
    if (grepl("/protein_id=", line)) {
      protein_id <- sub(".*/protein_id=\"([^\"]+)\".*", "\\1", line)
    }
    
    # If end of feature, add entry to list and reset variables
    if (grepl("^//", line) && !is.na(gene_id) && !is.na(locus_tag)) {
      annotation_list <- append(annotation_list, list(data.frame(GeneID = gene_id,
                                                                 locus_tag = locus_tag,
                                                                 transcript_id = transcript_id,
                                                                 protein_id = protein_id,
                                                                 stringsAsFactors = FALSE)))
      gene_id <- locus_tag <- transcript_id <- protein_id <- NA
    }
  }
  
  # Combine list into a data frame
  annotation_df <- bind_rows(annotation_list)
  
  # Join with the input data_frame on geneId
  merged_df <- data_frame %>%
    left_join(annotation_df, by = c("geneId" = "locus_tag"))
  
  # Save to CSV
  write.csv(merged_df, file = output_csv, row.names = FALSE)
  
  message("CSV file created successfully at: ", output_csv)
}

extract_annotations(anno.df_BbMSN2, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbMSN2_peaks_ID.csv")
extract_annotations(anno.df_BbNDT80, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbNDT80_peaks_ID.csv")
extract_annotations(anno.df_BbHCR1, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbHCR1_peaks_ID.csv")
extract_annotations(anno.df_BbFOXK2, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbFOXK2_peaks_ID.csv")
extract_annotations(anno.df_BbPacC, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbPacC_peaks_ID.csv")
extract_annotations(anno.df_BbOtf1, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbOft1_peaks_ID.csv")
extract_annotations(anno.df_BbCreA, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbCreA_peaks_ID.csv")
extract_annotations(anno.df_BbTBP, "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/genomic.gbff", "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbTBP_peaks_ID.csv")

#Filtrar intersección
library(purrr)

#CreA_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/CreA_proteinas.tsv"
BbMSN2_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbMSN2_proteinas.tsv"
BbNDT80_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbNDT80_proteinas.tsv"
BbPacC_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbPacC_proteinas.tsv"
BbFOXK2_proteinas_tsv <-  "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbFOXK2_proteinas.tsv"
BbOtf1_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbOtf1_proteinas.tsv"
BbHCR1_proteinas_tsv <- "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/BbHCR1_proteinas.tsv"

#df_CreA_proteinas <- read_tsv(CreA_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbMSN2_proteinas <- read_tsv(BbMSN2_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbNDT80_proteinas <- read_tsv(BbNDT80_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbPacC_proteinas <- read_tsv(BbPacC_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbFOXK2_proteinas <- read_tsv(BbFOXK2_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbOtf1_proteinas <- read_tsv(BbOtf1_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_BbHCR1_proteinas <- read_tsv(BbHCR1_proteinas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")

df_FT_proteinas_list <- list(df_BbMSN2_proteinas,df_BbNDT80_proteinas,df_BbPacC_proteinas,df_BbFOXK2_proteinas, df_BbOtf1_proteinas,df_BbHCR1_proteinas)
df_FT_proteinas_list_2 <- list(df_BbMSN2_proteinas,df_BbNDT80_proteinas,df_BbFOXK2_proteinas, df_BbOtf1_proteinas,df_BbHCR1_proteinas)
df_FT_proteinas_list_3 <- list(df_BbMSN2_proteinas,df_BbFOXK2_proteinas,df_BbHCR1_proteinas)

# Obtener intersección de locus_tag en todos los archivos
locus_comunes <- reduce(map(df_FT_proteinas_list, ~ .$locus_tag), intersect)
locus_comunes_2 <- reduce(map(df_FT_proteinas_list_2, ~ .$locus_tag), intersect)
locus_comunes_3 <- reduce(map(df_FT_proteinas_list_3, ~ .$locus_tag), intersect)

# Filtrar los dataframes para conservar solo los locus_tag en común
proteinas_filtradas <- map(df_FT_proteinas_list, ~ filter(.x, locus_tag %in% locus_comunes))
proteinas_filtradas_2 <- map(df_FT_proteinas_list_2, ~ filter(.x, locus_tag %in% locus_comunes_2))
proteinas_filtradas_3 <- map(df_FT_proteinas_list_3, ~ filter(.x, locus_tag %in% locus_comunes_3))

# Unir todos los datos en un solo dataframe
proteinas_filtradas_final <- bind_rows(proteinas_filtradas)
proteinas_filtradas_final_2 <- bind_rows(proteinas_filtradas_2)
proteinas_filtradas_final_3 <- bind_rows(proteinas_filtradas_3)

# Guardar en un nuevo archivo .tsv
write_tsv(proteinas_filtradas_final, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_coexpresadas.tsv")
write_tsv(proteinas_filtradas_final_2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_coexpresadas_sin_BbPacC.tsv")
write_tsv(proteinas_filtradas_final_3, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_expresadas_3.tsv")

Proteinas_coexpresadas_tsv <-  "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_coexpresadas.tsv"
Proteinas_coexpresadas_2_tsv <-  "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_coexpresadas_sin_BbPacC.tsv"
Proteinas_coexpresadas_3_tsv <-  "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_expresadas_3.tsv"

df_Proteinas_coexpresadas <- read_tsv(Proteinas_coexpresadas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_Proteinas_coexpresadas_2 <- read_tsv(Proteinas_coexpresadas_2_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_Proteinas_coexpresadas_3 <- read_tsv(Proteinas_coexpresadas_3_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
#Filtrar proteínas en .fasta de las que son predichas comoe expresión del FTs

# Extraer los protein_id del archivo .tsv
protein_BbMSN2_ids <- df_BbMSN2_proteinas$protein_id
protein_BbFOXK2_ids <- df_BbFOXK2_proteinas$protein_id
protein_BbNDT80_ids <- df_BbNDT80_proteinas$protein_id
protein_BbPacC_ids <- df_BbPacC_proteinas$protein_id
protein_BbOtf1_ids <- df_BbOtf1_proteinas$protein_id
protein_BbHCR1_ids <- df_BbHCR1_proteinas$protein_id
protein_coexpresadas_ids <- df_Proteinas_coexpresadas$protein_id
protein_coexpresadas_2_ids <- df_Proteinas_coexpresadas_2$protein_id
protein_coexpresadas_3_ids <- df_Proteinas_coexpresadas_3$protein_id

# Paso 2: Leer el archivo FASTA con las proteínas
proteinas_fasta <- readAAStringSet("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/protein.faa")
nombres_proteinas_fasta <- names(proteinas_fasta)

# Paso 3: Filtrar las proteínas que están en el .tsv usando protein_id
# Convertimos las cabeceras del FASTA (que contienen los protein_id) en un vector
fasta_ids <- sapply(strsplit(names(proteinas_fasta), " "), `[`, 1)

# Filtramos las proteínas cuyos protein_id están en el archivo .tsv
proteinas_filtradas_BbFOXK2 <- proteinas_fasta[fasta_ids %in% protein_BbFOXK2_ids]
proteinas_filtradas_BbNDT80 <- proteinas_fasta[fasta_ids %in% protein_BbNDT80_ids]
proteinas_filtradas_BbPacC <- proteinas_fasta[fasta_ids %in% protein_BbPacC_ids]
proteinas_filtradas_BbMSN2 <- proteinas_fasta[fasta_ids %in% protein_BbMSN2_ids]
proteinas_filtradas_BbOtf1 <- proteinas_fasta[fasta_ids %in% protein_BbOtf1_ids]
proteinas_filtradas_BbHCR1 <- proteinas_fasta[fasta_ids %in% protein_BbHCR1_ids]
proteinas_filtradas_coexpresadas <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_ids]
proteinas_filtradas_coexpresadas_2 <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_2_ids]
proteinas_filtradas_coexpresadas_3 <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_3_ids]

# Paso 4: Guardar el resultado en un nuevo archivo FASTA
writeXStringSet(proteinas_filtradas_BbFOXK2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbFOXK2.fasta")
writeXStringSet(proteinas_filtradas_BbNDT80, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbNDT80.fasta")
writeXStringSet(proteinas_filtradas_BbPacC, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbPacC.fasta")
writeXStringSet(proteinas_filtradas_BbMSN2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbMSN2.fasta")
writeXStringSet(proteinas_filtradas_BbOtf1, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbOtf1.fasta")
writeXStringSet(proteinas_filtradas_BbHCR1, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbHCR1.fasta")
writeXStringSet(proteinas_filtradas_coexpresadas, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_coexpresadas.fasta")
writeXStringSet(proteinas_filtradas_coexpresadas_2, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_coexpresadas_2.fasta")
writeXStringSet(proteinas_filtradas_coexpresadas_3, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_coexpresadas_3.fasta")

#Filtrar por Signal_P
library(readr)
library(dplyr)

# Leer el archivo ignorando las líneas que comienzan con #
SignalP_data_BbFOXK2 <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbFOXK2.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_BbHCR1 <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbHCR1.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_BbPacC <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbPacC.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_BbMSN2 <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbMSN2.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_BbNDT80 <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbNDT80.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_BbOtf1 <- read_tsv("/home/Nephelim/Descargas/Signal_P_BbOtf1.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_Coexpresion_2 <- read_tsv("/home/Nephelim/Descargas/Signal_P_Coexpresion_2.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_Coexpresion_3 <- read_tsv("/home/Nephelim/Descargas/Signal_P_Coexpresion_3.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)

colnames(SignalP_data_BbFOXK2) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_BbHCR1) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_BbPacC) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_BbMSN2) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_BbNDT80) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_BbOtf1) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_Coexpresion_2) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")
colnames(SignalP_data_Coexpresion_3) <- c("ID", "Prediction", "SP_Sec_SPI", "OTHER", "CS_Position")

SignalP_data_BbFOXK2 <- SignalP_data_BbFOXK2 %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbFOXK2 <- SignalP_data_BbFOXK2 %>% filter(!is.na(CS_Position))


SignalP_data_BbHCR1 <- SignalP_data_BbHCR1 %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbHCR1 <- SignalP_data_BbHCR1 %>% filter(!is.na(CS_Position))

SignalP_data_BbPacC <- SignalP_data_BbPacC %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbPacC <- SignalP_data_BbPacC %>% filter(!is.na(CS_Position))

SignalP_data_BbMSN2 <- SignalP_data_BbMSN2 %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbMSN2 <- SignalP_data_BbMSN2 %>% filter(!is.na(CS_Position))

SignalP_data_BbNDT80 <- SignalP_data_BbNDT80 %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbNDT80 <- SignalP_data_BbNDT80 %>% filter(!is.na(CS_Position))

SignalP_data_BbOtf1 <- SignalP_data_BbOtf1 %>% filter(!is.na(ID) & ID != "")
SignalP_data_BbOtf1 <- SignalP_data_BbOtf1 %>% filter(!is.na(CS_Position))

SignalP_data_Coexpresion_2 <- SignalP_data_Coexpresion_2 %>% filter(!is.na(ID) & ID != "")
SignalP_data_Coexpresion_2 <- SignalP_data_Coexpresion_2 %>% filter(!is.na(CS_Position))

SignalP_data_Coexpresion_3 <- SignalP_data_Coexpresion_3 %>% filter(!is.na(ID) & ID != "")
SignalP_data_Coexpresion_3 <- SignalP_data_Coexpresion_3 %>% filter(!is.na(CS_Position))

# Ver los primeros registros
head(SignalP_data_BbFOXK2)
head(SignalP_data_BbHCR1)
head(SignalP_data_BbPacC)
head(SignalP_data_BbMSN2)
head(SignalP_data_BbNDT80)
head(SignalP_data_BbOtf1)
head(SignalP_data_Coexpresion_2)
head(SignalP_data_Coexpresion_3)

# Pasar a data.frame
SignalP_data_BbFOXK2_df <- as.data.frame(SignalP_data_BbFOXK2)
SignalP_data_BbHCR1_df <- as.data.frame(SignalP_data_BbHCR1)
SignalP_data_BbPacC_df <- as.data.frame(SignalP_data_BbPacC)
SignalP_data_BbMSN2_df <- as.data.frame(SignalP_data_BbMSN2)
SignalP_data_BbNDT80_df <- as.data.frame(SignalP_data_BbNDT80)
SignalP_data_BbOtf1_df <- as.data.frame(SignalP_data_BbOtf1)
SignalP_data_Coexpresion_2_df <- as.data.frame(SignalP_data_Coexpresion_2)
SignalP_data_Coexpresion_3_df <- as.data.frame(SignalP_data_Coexpresion_3)

# Obtener los protein_id
proteinas_id_BbFOXK2_Signal_P <- SignalP_data_BbFOXK2_df$ID
proteinas_id_BbHCR1_Signal_P <- SignalP_data_BbHCR1_df$ID
proteinas_id_BbPacC_Signal_P <- SignalP_data_BbPacC_df$ID
proteinas_id_BbMSN2_Signal_P <- SignalP_data_BbMSN2_df$ID
proteinas_id_BbNDT80_Signal_P <- SignalP_data_BbNDT80_df$ID
proteinas_id_BbOtf1_Signal_P <- SignalP_data_BbOtf1_df$ID
proteinas_id_Coexpresion_2_Signal_P <- SignalP_data_Coexpresion_2_df$ID
proteinas_id_Coexpresion_3_Signal_P <- SignalP_data_Coexpresion_3_df$ID

# Filtrar las proteinas cuyosprotein_id están en el archivo .tsv

proteinas_filtradas_BbFOXK2_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbFOXK2_Signal_P]
proteinas_filtradas_BbHCR1_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbHCR1_Signal_P]
proteinas_filtradas_BbPacC_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbPacC_Signal_P]
proteinas_filtradas_BbMSN2_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbMSN2_Signal_P]
proteinas_filtradas_BbNDT80_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbNDT80_Signal_P]
proteinas_filtradas_BbOtf1_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_BbOtf1_Signal_P]
proteinas_filtradas_Coexpresion_2_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_Coexpresion_2_Signal_P]
proteinas_filtradas_Coexpresion_3_Signal_P <- proteinas_fasta[fasta_ids %in% proteinas_id_Coexpresion_3_Signal_P]

#Guardar el resultado en un archivo .FASTA

writeXStringSet(proteinas_filtradas_BbFOXK2_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbFOXK2_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_BbHCR1_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbHCR1_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_BbPacC_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbPacC_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_BbMSN2_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbMSN2_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_BbNDT80_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbNDT80_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_BbOtf1_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_BbOtf1_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_Coexpresion_2_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_Coexpresion_2_Signal_P.fasta")
writeXStringSet(proteinas_filtradas_Coexpresion_3_Signal_P, "/home/Nephelim/Descargas/Factores de transcripcion/Transcribed proteins/proteinas_filtradas_Coexpresion_3_Signal_P.fasta")

#Análisis Enriquecimiento - GO


# Leer archivo GAF <- Contiene los GO annotations
gaf_file <- "/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/GCF_000280675.1_ASM28067v1_gene_ontology.gaf"
gaf_data <- read.delim(gaf_file, header = FALSE, comment.char = "!", sep = "\t", stringsAsFactors = FALSE)

# Extraer columnas necesarias: ID del gen (DB_Object_ID) y términos GO (GO_ID)
colnames(gaf_data) <- c("DB", "GeneID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", "With_Or_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By")

# Verificar las primeras filas
head(gaf_data)

# Lista de genes anotados en los picos
# Suponiendo que tu archivo annotated-peaks tiene una columna `gene_id` con los IDs de NCBI
BbCreA_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbCreA.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbCreA_gene_list <- unique(BbCreA_genes_proteins_peeaks[["GeneID"]])
BbCreA_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbCreA_gene_list))

BbMSN2_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbMSN2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbMSN2_gene_list <- unique(BbMSN2_genes_proteins_peeaks[["GeneID"]])
BbMSN2_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbMSN2_gene_list))

BbFOXK2_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbFOXK2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbFOXK2_gene_list <- unique(BbFOXK2_genes_proteins_peeaks[["GeneID"]])
BbFOXK2_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbFOXK2_gene_list))

BbHCR1_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbHCR1.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbHCR1_gene_list <- unique(BbHCR1_genes_proteins_peeaks[["GeneID"]])
BbHCR1_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbHCR1_gene_list))

BbNDT80_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbNDT80.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbNDT80_gene_list <- unique(BbNDT80_genes_proteins_peeaks[["GeneID"]])
BbNDT80_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbNDT80_gene_list))

BbOtf1_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbOtf1.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbOtf1_gene_list <- unique(BbOtf1_genes_proteins_peeaks[["GeneID"]])
BbOtf1_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbOtf1_gene_list))

BbPacC_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbPacC.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbPacC_gene_list <- unique(BbPacC_genes_proteins_peeaks[["GeneID"]])
BbPacC_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbPacC_gene_list))

BbTBP_genes_proteins_peeaks <- read.table("/home/Nephelim/Descargas/Factores de transcripcion/Annotated Peaks/anno.df_BbTBP.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
BbTBP_gene_list <- unique(BbTBP_genes_proteins_peeaks[["GeneID"]])
BbTBP_gene_selection <- factor(as.integer(gaf_data$GeneID %in% BbTBP_gene_list))

names(BbMSN2_gene_selection) <- gaf_data$GeneID
names(BbCreA_gene_selection) <- gaf_data$GeneID
names(BbFOXK2_gene_selection) <- gaf_data$GeneID
names(BbHCR1_gene_selection) <- gaf_data$GeneID
names(BbNDT80_gene_selection) <- gaf_data$GeneID
names(BbPacC_gene_selection) <- gaf_data$GeneID
names(BbOtf1_gene_selection) <- gaf_data$GeneID
names(BbTBP_gene_selection) <- gaf_data$GeneID

# Crear un objeto topGOdata para análisis de enriquecimiento
GOdata_BbCreA_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbCreA_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbCreA_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbCreA_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbCreA_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbCreA_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))

GOdata_BbMSN2_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbMSN2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbMSN2_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbMSN2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbMSN2_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbMSN2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))

GOdata_BbFOXK2_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbFOXK2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbFOXK2_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbFOXK2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbFOXK2_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbFOXK2_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))


GOdata_BbHCR1_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbHCR1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbHCR1_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbHCR1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbHCR1_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbHCR1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))


GOdata_BbNDT80_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbNDT80_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbNDT80_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbNDT80_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbNDT80_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbNDT80_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))

GOdata_BbOtf1_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbOtf1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbOtf1_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbOtf1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbOtf1_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbOtf1_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))

GOdata_BbPacC_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbPacC_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbPacC_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbPacC_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbPacC_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbPacC_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))


GOdata_BbTBP_BP <- new("topGOdata",
              ontology = "BP",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbTBP_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbTBP_MF <- new("topGOdata",
              ontology = "MF",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbTBP_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))
GOdata_BbTBP_CC <- new("topGOdata",
              ontology = "CC",  # O usa "MF" o "CC" si prefieres otros aspectos
              allGenes = BbTBP_gene_selection,
              nodeSize = 10, # nodeSize filtra términos GO con pocos genes
              annot = annFUN.gene2GO,
              gene2GO = split(gaf_data$GO_ID, gaf_data$GeneID))


# Realizar el análisis usando el método 'classic' o 'elim' y el test estadístico Fisher
result_fisher_GO_BbCreA_BP <- runTest(GOdata_BbCreA_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbCreA_MF <- runTest(GOdata_BbCreA_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbCreA_CC <- runTest(GOdata_BbCreA_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbCreA_BP <- GenTable(GOdata_BbCreA_BP,
                        classicFisher = result_fisher_GO_BbCreA_BP,
                        orderBy = "classicFisher", 
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbCreA_MF <- GenTable(GOdata_BbCreA_MF,
                        classicFisher = result_fisher_GO_BbCreA_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbCreA_CC <- GenTable(GOdata_BbCreA_CC,
                        classicFisher = result_fisher_GO_BbCreA_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbMSN2_BP <- runTest(GOdata_BbMSN2_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbMSN2_MF <- runTest(GOdata_BbMSN2_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbMSN2_CC <- runTest(GOdata_BbMSN2_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbMSN2_BP <- GenTable(GOdata_BbMSN2_BP,
                        classicFisher = result_fisher_GO_BbMSN2_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbMSN2_MF <- GenTable(GOdata_BbMSN2_MF,
                        classicFisher = result_fisher_GO_BbMSN2_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbMSN2_CC <- GenTable(GOdata_BbMSN2_CC,
                        classicFisher = result_fisher_GO_BbMSN2_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbFOXK2_BP <- runTest(GOdata_BbFOXK2_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbFOXK2_MF <- runTest(GOdata_BbFOXK2_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbFOXK2_CC <- runTest(GOdata_BbFOXK2_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbFOXK2_BP <- GenTable(GOdata_BbFOXK2_BP,
                        classicFisher = result_fisher_GO_BbFOXK2_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbFOXK2_MF <- GenTable(GOdata_BbFOXK2_MF,
                        classicFisher = result_fisher_GO_BbFOXK2_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbFOXK2_CC <- GenTable(GOdata_BbFOXK2_CC,
                        classicFisher = result_fisher_GO_BbFOXK2_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbHCR1_BP <- runTest(GOdata_BbHCR1_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbHCR1_MF <- runTest(GOdata_BbHCR1_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbHCR1_CC <- runTest(GOdata_BbHCR1_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbHCR1_BP <- GenTable(GOdata_BbHCR1_BP,
                        classicFisher = result_fisher_GO_BbHCR1_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbHCR1_MF <- GenTable(GOdata_BbHCR1_MF,
                        classicFisher = result_fisher_GO_BbHCR1_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbHCR1_CC <- GenTable(GOdata_BbHCR1_CC,
                        classicFisher = result_fisher_GO_BbHCR1_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbNDT80_BP <- runTest(GOdata_BbNDT80_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbNDT80_MF <- runTest(GOdata_BbNDT80_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbNDT80_CC <- runTest(GOdata_BbNDT80_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbNDT80_BP <- GenTable(GOdata_BbNDT80_BP,
                        classicFisher = result_fisher_GO_BbNDT80_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbNDT80_MF <- GenTable(GOdata_BbNDT80_MF,
                        classicFisher = result_fisher_GO_BbNDT80_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbNDT80_CC <- GenTable(GOdata_BbNDT80_CC,
                        classicFisher = result_fisher_GO_BbNDT80_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbOtf1_BP <- runTest(GOdata_BbOtf1_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbOtf1_MF <- runTest(GOdata_BbOtf1_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbOtf1_CC <- runTest(GOdata_BbOtf1_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbOtf1_BP <- GenTable(GOdata_BbOtf1_BP,
                        classicFisher = result_fisher_GO_BbOtf1_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbOtf1_MF <- GenTable(GOdata_BbOtf1_MF,
                        classicFisher = result_fisher_GO_BbOtf1_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbOtf1_CC <- GenTable(GOdata_BbOtf1_CC,
                        classicFisher = result_fisher_GO_BbOtf1_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados

result_fisher_GO_BbPacC_BP <- runTest(GOdata_BbPacC_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbPacC_MF <- runTest(GOdata_BbPacC_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbPacC_CC <- runTest(GOdata_BbPacC_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbPacC_BP <- GenTable(GOdata_BbPacC_BP,
                        classicFisher = result_fisher_GO_BbPacC_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbPacC_MF <- GenTable(GOdata_BbPacC_MF,
                        classicFisher = result_fisher_GO_BbPacC_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbPacC_CC <- GenTable(GOdata_BbPacC_CC,
                        classicFisher = result_fisher_GO_BbPacC_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados


result_fisher_GO_BbTBP_BP <- runTest(GOdata_BbTBP_BP, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbTBP_MF <- runTest(GOdata_BbTBP_MF, algorithm = "classic", statistic = "fisher")
result_fisher_GO_BbTBP_CC <- runTest(GOdata_BbTBP_CC, algorithm = "classic", statistic = "fisher")

top_results_GO_BbTBP_BP <- GenTable(GOdata_BbTBP_BP,
                        classicFisher = result_fisher_GO_BbTBP_BP,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbTBP_MF <- GenTable(GOdata_BbTBP_MF,
                        classicFisher = result_fisher_GO_BbTBP_MF,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados
top_results_GO_BbTBP_CC <- GenTable(GOdata_BbTBP_CC,
                        classicFisher = result_fisher_GO_BbTBP_CC,
                        orderBy = "classicFisher",
                        topNodes = 10)  # Ajusta topNodes según el número de resultados deseados


# Crear un gráfico de barras con los términos GO enriquecidos
ggplot(top_results_GO_BbCreA_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "cadetblue3") +
  labs(title = "Principales términos GO de BbCreA (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_CreA_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbCreA_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "cadetblue3") +
  labs(title = "Principales términos GO de BbCreA (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_CreA_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbCreA_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "cadetblue3") +
  labs(title = "Principales términos GO de BbCreA (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_CreA_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbMSN2_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1A1A") +
  labs(title = "Principales términos GO de BbMSN2 (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbMSN2_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbMSN2_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1A1A") +
  labs(title = "Principales términos GO de BbMSN2 (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbMSN2_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbMSN2_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1A1A") +
  labs(title = "Principales términos GO de BbMSN2 (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbMSN2_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbFOXK2_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#008B00") +
  labs(title = "Principales términos GO de BbFOXK2 (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbFOXK2_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbFOXK2_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#008B00") +
  labs(title = "Principales términos GO de BbFOXK2 (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbFOXK2_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbFOXK2_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#008B00") +
  labs(title = "Principales términos GO de BbFOXK2(Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbFOXK2_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbPacC_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1C62") +
  labs(title = "Principales términos GO de BbPacC (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbPacC_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbPacC_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1C62") +
  labs(title = "Principales términos GO de BbPacC (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbPacC_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbPacC_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#8B1C62") +
  labs(title = "Principales términos GO de BbPacC (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbPacC_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)


ggplot(top_results_GO_BbHCR1_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#5D478B") +
  labs(title = "Principales términos GO de BbHCR1 (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbHCR1_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbHCR1_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#5D478B") +
  labs(title = "Principales términos GO de BbHCR1 (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbHCR1_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbHCR1_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#5D478B") +
  labs(title = "Principales términos GO de BbHCR1 (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbHCR1_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)


ggplot(top_results_GO_BbNDT80_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "palevioletred4") +
  labs(title = "Principales términos GO de BbNDT80 (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbNDT80_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbNDT80_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "palevioletred4") +
  labs(title = "Principales términos GO de BbNDT80 (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbNDT80_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbNDT80_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "palevioletred4") +
  labs(title = "Principales términos GO de BbNDT80 (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbNDT80_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)


ggplot(top_results_GO_BbOtf1_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#27408B") +
  labs(title = "Principales términos GO de BbOtf1 (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbOtf1_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbOtf1_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#27408B") +
  labs(title = "Principales términos GO de BbOtf1 (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbOtf1_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbOtf1_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "#27408B") +
  labs(title = "Principales términos GO de BbOtf1 (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbOtf1_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)


ggplot(top_results_GO_BbTBP_BP, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "tan3") +
  labs(title = "Principales términos GO de BbTBP (Procesos Biológicos)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_minimal(base_size = 14) +  # Aumenta el tamaño general del texto
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) + # Añade más espacio a la izquierda
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbTBP_BP_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbTBP_MF, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "tan3") +
  labs(title = "Principales términos GO de BbTBP (Función molecular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbTBP_MF_barplot_fixed.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbTBP_CC, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity", fill = "tan3") +
  labs(title = "Principales términos GO de BbTBP (Compartimiento celular)",
       x = "Término GO", y = "-log10(p-value)") +
  coord_flip() +
  theme_classic(base_size = 14) +  # Utiliza theme_classic para un aspecto más limpio
  theme(plot.title = element_text(hjust = 0.15, size = 36, face = "bold", family = "Helvetica"),
        strip.text = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 24, angle = 45, face = "bold", hjust = 1),  # Rotar etiquetas
        axis.text.y = element_text(family="Helvetica", size = 24),
        axis.title.x = element_text(family = "Helvetica", size = 28, face = "bold"),
        axis.title.y = element_text(family="Helvetica", size = 28, face = "bold"),
        plot.margin = margin(1, 1, 1, 3, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
ggsave("GO_BbTBP_CC_barplot_fixed.png", width = 24, height = 16, dpi = 300)



# Asegurarse de tener la columna del recuento de genes en los términos GO enriquecidos
top_results_GO_BbCreA_BP$Significant <- as.numeric(top_results_GO_BbCreA_BP$Significant)  # Convertir el número de genes asociados a numérico
top_results_GO_BbCreA_MF$Significant <- as.numeric(top_results_GO_BbCreA_MF$Significant)
top_results_GO_BbCreA_CC$Significant <- as.numeric(top_results_GO_BbCreA_CC$Significant)

ggplot(top_results_GO_BbCreA_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbCreA - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbCreA_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbCreA_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbCreA - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbCreA_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)


ggplot(top_results_GO_BbCreA_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbCreA - Compartimiento Celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbCreA_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbMSN2_BP$Significant <- as.numeric(top_results_GO_BbMSN2_BP$Significant)
top_results_GO_BbMSN2_MF$Significant <- as.numeric(top_results_GO_BbMSN2_MF$Significant)
top_results_GO_BbMSN2_CC$Significant <- as.numeric(top_results_GO_BbMSN2_CC$Significant)

ggplot(top_results_GO_BbMSN2_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbMSN2 - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbMSN2_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbMSN2_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbMSN2 - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbMSN2_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbMSN2_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbMSN2 - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbMSN2_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbFOXK2_BP$Significant <- as.numeric(top_results_GO_BbFOXK2_BP$Significant)
top_results_GO_BbFOXK2_MF$Significant <- as.numeric(top_results_GO_BbFOXK2_MF$Significant)
top_results_GO_BbFOXK2_CC$Significant <- as.numeric(top_results_GO_BbFOXK2_CC$Significant)

ggplot(top_results_GO_BbFOXK2_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbFOXK2 - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbFOXK2_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbFOXK2_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbFOXK2 - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbFOXK2_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbFOXK2_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbFOXK2 - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbFOXK2_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbNDT80_BP$Significant <- as.numeric(top_results_GO_BbNDT80_BP$Significant)
top_results_GO_BbNDT80_MF$Significant <- as.numeric(top_results_GO_BbNDT80_MF$Significant)
top_results_GO_BbNDT80_CC$Significant <- as.numeric(top_results_GO_BbNDT80_CC$Significant)

ggplot(top_results_GO_BbNDT80_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbNDT80 - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbNDT80_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbNDT80_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbNDT80 - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbNDT80_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbNDT80_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbNDT80 - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbNDT80_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbHCR1_BP$Significant <- as.numeric(top_results_GO_BbHCR1_BP$Significant)
top_results_GO_BbHCR1_MF$Significant <- as.numeric(top_results_GO_BbHCR1_MF$Significant)
top_results_GO_BbHCR1_CC$Significant <- as.numeric(top_results_GO_BbHCR1_CC$Significant)

ggplot(top_results_GO_BbHCR1_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbHCR1 - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbHCR1_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbHCR1_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbHCR1 - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbHCR1_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbHCR1_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbHCR1 - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbHCR1_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbPacC_BP$Significant <- as.numeric(top_results_GO_BbPacC_BP$Significant)
top_results_GO_BbPacC_MF$Significant <- as.numeric(top_results_GO_BbPacC_MF$Significant)
top_results_GO_BbPacC_CC$Significant <- as.numeric(top_results_GO_BbPacC_CC$Significant)

ggplot(top_results_GO_BbPacC_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbPacC - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbPacC_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbPacC_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbPacC - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbPacC_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbPacC_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbPacC - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbPacC_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbOtf1_BP$Significant <- as.numeric(top_results_GO_BbOtf1_BP$Significant)
top_results_GO_BbOtf1_MF$Significant <- as.numeric(top_results_GO_BbOtf1_MF$Significant)
top_results_GO_BbOtf1_CC$Significant <- as.numeric(top_results_GO_BbOtf1_CC$Significant)

ggplot(top_results_GO_BbOtf1_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbOtf1 - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbOtf1_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbOtf1_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbOtf1 - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbOtf1_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbOtf1_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbOtf1 - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbOtf1_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)

top_results_GO_BbTBP_BP$Significant <- as.numeric(top_results_GO_BbTBP_BP$Significant)
top_results_GO_BbTBP_MF$Significant <- as.numeric(top_results_GO_BbTBP_MF$Significant)
top_results_GO_BbTBP_CC$Significant <- as.numeric(top_results_GO_BbTBP_CC$Significant)

ggplot(top_results_GO_BbTBP_BP, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbTBP - Procesos Biológicos") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbTBP_BP_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbTBP_MF, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
 geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbTBP - Función Molecular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbTBP_MF_dotplot_significant.png", width = 24, height = 16, dpi = 300)

ggplot(top_results_GO_BbOtf1_CC, aes(x = Significant, y = reorder(Term, -as.numeric(classicFisher)), color = as.numeric(classicFisher))) +
  geom_point(aes(size = Significant), stroke = 1.5) +
  scale_color_gradient(low = "darkorchid3", high = "coral") +
  scale_size_continuous(range = c(4, 20)) +
  labs(x = "Número de genes", y = "Términos GO", color = "p-Valor (Fisher)", size = "Número de genes") +
  ggtitle("Gráfico Dispersión de Términos GO enriquecidos - BbTBP - Compartimiento celular") +
  theme_bw() +  # Aumenta el tamaño general del texto
  theme(
      plot.title = element_text(hjust = 0.75, size = 37, face = "bold", family = "Helvetica", color = "#4B4B4D"),
      axis.text.x = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      axis.text.y = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold", color = "#4B4B4D"),
      axis.title.y = element_text(family = "Helvetica", size = 30, face = "bold", color = "#4B4B4D"),
      legend.title = element_text(family = "Helvetica", size = 21, face = "bold", color = "#4B4B4D"),
      legend.text = element_text(family = "Helvetica", size = 19, face = "bold", color = "#4B4B4D"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Agregar borde al gráfico
      plot.margin = margin(10, 30, 10, 30),  # Ajustar márgenes para dar espacio
      axis.ticks.length = unit(0.3, "cm"))  # Aumentar longitud de los ticks
ggsave("GO_BbTBP_CC_dotplot_significant.png", width = 24, height = 16, dpi = 300)
