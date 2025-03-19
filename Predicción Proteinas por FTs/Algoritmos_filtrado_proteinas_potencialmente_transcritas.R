#Si no se tienen las librerías:

#install.packages("purrr")
#install.packages("dplyr")
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("ggplot2")

library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)


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
extract_geneid_locus_tag(ruta_archivo_gbff, ruta_archivo_gene_mapping_tsv_1)

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

extract_gene_info_from_gtf(ruta_archivo_gtf, ruta_archivo_gene_mapping_tsv_2)

gene_mapping <- read.table(ruta_archivo_gene_mapping_tsv_1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_mapping_locusTagOnly <- read.table(ruta_archivo_gene_mapping_tsv_2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

head(gene_mapping)
head(gene_mapping_locusTagOnly)

########

#Cargar los Data_frame de los elementos transcritos por los FTs:

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

#Correcciones a los DataFrame
anno_dfs_2 <- list()

for (factor in factores_transcripcion){
  anno_dfs_2[[factor]] <- merge(anno_dfs[[factor]], gene_mapping, by.x = "geneId", by.y = "locus_tag", all.x = TRUE)
  anno_dfs_2[[factor]] <- anno_dfs_2[[factor]]%>% rename(locus_tag = geneId)
}

ruta_anotaciones_fimo <- ""

for (factor in factores_transcripcion){
  write.table(anno_dfs_2[[factor]], paste(ruta_anotaciones_fimo, "anno.df_", factor, ".tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(anno_dfs[[factor]], paste(ruta_anotaciones_fimo, "anno.df_", factor, "_origin.tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#Comprobar lo guardado
anno_dfs <- list()

for (factor in factores_transcripcion){
  anno_dfs[[factor]] <- read.table(paste(ruta_anotaciones_fimo, "anno.df_", factor, ".tsv", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

# Extraer los genes proteicos
genesIDs <- list()

for (factor in factores_transcripcion){
  genesIDs[[factor]] <- anno_dfs[[factor]] %>% select(locus_tag, GeneID, transcript_id, protein_id, product)
  genesIDs[[factor]] <- na.omit(genesIDs[[factor]])
  genesIDs[[factor]] <- unique(genesIDs[[factor]])
}


# Crear un gráfico de barras para comparar los factores
gene_count <- c()

for (factor in factores_transcripcion){
  gene_count <- c(gene_count, length(genesIDs[[factor]]$GeneID))
}

# Dataframe para el gráfico
df_genes_FTs <- data.frame(Factor = factor(factores_transcripcion, levels = factores_transcripcion), Count = gene_count)

# Gráfico de barras mejorado
ggplot(df_genes_FTs, aes(x = Factor, y = Count, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "plasma") +  # Paleta con más colores
  labs(x = "Factor de Transcripción", y = "Número de genes transcritos",
       title = "Comparación de Cantidad de Genes entre Factores de Transcripción") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 25, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 25),
        axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 22, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank()) +
  coord_flip()  # Cambia la orientación del gráfico para mostrar mejor los



#Generar un .tsv de cada Factor de transcripción que contenga las proteínas las cuales induce su transcripción.
ruta_proteinas_transcriptas <- ""

for (factor in factores_transcripcion){
  write.table(genesIDs[[factor]], paste(ruta_proteinas_transcriptas, factor, "_proteinas.tsv", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#Otra función para generar un .csv que tenga los annotation peaks con los valores de annotacion protein_id, transcript_id, locus_tag, etc.
extract_annotations <- function(data_frame, gbff_path, output_csv) {
  
  # Leer el archivo .gbff 
  gbff_content <- readLines(gbff_path)
  
  # "parsear" las anotaciones del archivo .gbff 
  annotation_list <- list()
  gene_id <- locus_tag <- transcript_id <- protein_id <- NA
  
  for (line in gbff_content) {
    
    # Chequear GeneID
    if (grepl("/db_xref=\"GeneID:", line)) {
      gene_id <- sub(".*/db_xref=\"GeneID:(\\d+)\".*", "\\1", line)
    }
    
    # Chequear locus_tag
    if (grepl("/locus_tag=", line)) {
      locus_tag <- sub(".*/locus_tag=\"([^\"]+)\".*", "\\1", line)
    }
    
    # Chequear transcript_id
    if (grepl("/transcript_id=", line)) {
      transcript_id <- sub(".*/transcript_id=\"([^\"]+)\".*", "\\1", line)
    }
    
    # Chequear protein_id si existe
    if (grepl("/protein_id=", line)) {
      protein_id <- sub(".*/protein_id=\"([^\"]+)\".*", "\\1", line)
    }
    
    # Si se llega al fnal de la características, añadir cada entrada a la lista y resetear las variables
    if (grepl("^//", line) && !is.na(gene_id) && !is.na(locus_tag)) {
      annotation_list <- append(annotation_list, list(data.frame(GeneID = gene_id,
                                                                 locus_tag = locus_tag,
                                                                 transcript_id = transcript_id,
                                                                 protein_id = protein_id,
                                                                 stringsAsFactors = FALSE)))
      gene_id <- locus_tag <- transcript_id <- protein_id <- NA
    }
  }
  
  # Combinar lista en un data_frame
  annotation_df <- bind_rows(annotation_list)
  
  # Unir el data_farame por locus_tag
  merged_df <- data_frame %>%
    left_join(annotation_df, by = c("geneId" = "locus_tag"))
  
  # Guardar en CSV
  write.csv(merged_df, file = output_csv, row.names = FALSE)
  
  message("CSV file created successfully at: ", output_csv)
}

ruta_genomic_gbff <- ""

for(factor in factores_transcripcion){
  extract_annotations(anno_dfs[[factor]], ruta_genomic_gbff, pase(ruta_proteinas_transcriptas, factor, "_peaks_ID.csv", sep = ""))
}

FTs_proteinas_tsv <- list()

for(factor in factores_transcripcion){
  FTs_proteinas_tsv[[factor]] <- paste(ruta_proteinas_transcriptas, factor, "_proteinas.tsv", sep = "")
}

df_proteinas_tsv <- list()

for(factor in factores_transcripcion){
  df_proteinas_tsv[[factor]] <- read_tsv(FTs_proteinas_tsv[[factor]], show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
}

factores_seleccionados <- c("BbAP1", "BbCreA", "BbCrz1", "BbFlbC", "BbHcr1", "BbKlf", "BbOps3", "BbOsrR1", "BbOsrR3", "BbOtf1",  "BbRei1", "BbRep1") #Segun BP - GO comparten procesos
factores_seleccionados_2_1 <- c("BbCtf1alpha","BbHox2", "BbSmr1") #Segun BP - GO comparten procesos
factores_seleccionados_2_2 <- c("BbMcm1", "BbTcp1") #Segun BP - GO comparten procesos

df_proteinas_coexpresadas_2_tsv <- list()

for(factor in factores_seleccionados){
  df_proteinas_coexpresadas_2_tsv[[factor]] <- read_tsv(FTs_proteinas_tsv[[factor]], show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
}

df_proteinas_coexpresadas_3_1_tsv <- list()

for(factor in factores_seleccionados_2_1){
  df_proteinas_coexpresadas_3_1_tsv[[factor]] <- read_tsv(FTs_proteinas_tsv[[factor]], show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
}

df_proteinas_coexpresadas_3_2_tsv <- list()

for(factor in factores_seleccionados_2_2){
  df_proteinas_coexpresadas_3_2_tsv[[factor]] <- read_tsv(FTs_proteinas_tsv[[factor]], show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
}

# Obtener intersección de locus_tag en todos los archivos
locus_comunes <- reduce(map(df_proteinas_tsv, ~ .$locus_tag), intersect)
locus_comunes_2 <- reduce(map(df_proteinas_coexpresadas_2_tsv, ~ .$locus_tag), intersect)
locus_comunes_3_1 <- reduce(map(df_proteinas_coexpresadas_3_1_tsv, ~ .$locus_tag), intersect)
locus_comunes_3_2 <- reduce(map(df_proteinas_coexpresadas_3_2_tsv, ~ .$locus_tag), intersect)

# Filtrar los dataframes para conservar solo los locus_tag en común
proteinas_filtradas <- map(df_proteinas_tsv, ~ filter(.x, locus_tag %in% locus_comunes))
proteinas_filtradas_2 <- map(df_proteinas_coexpresadas_2_tsv, ~ filter(.x, locus_tag %in% locus_comunes_2))
proteinas_filtradas_3_1 <- map(df_proteinas_coexpresadas_3_1_tsv, ~ filter(.x, locus_tag %in% locus_comunes_3_1))
proteinas_filtradas_3_2 <- map(df_proteinas_coexpresadas_3_2_tsv, ~ filter(.x, locus_tag %in% locus_comunes_3_2))

# Unir todos los datos en un solo dataframe
proteinas_filtradas_final <- bind_rows(proteinas_filtradas)
proteinas_filtradas_final <- unique(proteinas_filtradas_final)
proteinas_filtradas_final_2 <- bind_rows(proteinas_filtradas_2)
proteinas_filtradas_final_2 <- unique(proteinas_filtradas_final_2)
proteinas_filtradas_final_3_1 <- bind_rows(proteinas_filtradas_3_1)
proteinas_filtradas_final_3_1 <- unique(proteinas_filtradas_final_3_1)
proteinas_filtradas_final_3_2 <- bind_rows(proteinas_filtradas_3_2)
proteinas_filtradas_final_3_2 <- unique(proteinas_filtradas_final_3_2)

# Guardar en un nuevo archivo .tsv
write_tsv(proteinas_filtradas_final, "/home/Nephelim/Descargas/Proteinas trasncriptas/proteinas_coexpresadas.tsv")
write_tsv(proteinas_filtradas_final_2, "/home/Nephelim/Descargas/Proteinas trasncriptas/proteinas_coexpresadas_2.tsv")
write_tsv(proteinas_filtradas_final_3_1, "/home/Nephelim/Descargas/Proteinas trasncriptas/proteinas_coexpresadas_3_1.tsv")
write_tsv(proteinas_filtradas_final_3_2, "/home/Nephelim/Descargas/Proteinas trasncriptas/proteinas_coexpresadas_3_2.tsv")

Proteinas_coexpresadas_tsv <-  paste(ruta_proteinas_transcriptas, "proteinas_coexpresadas.tsv", sep = "")
Proteinas_coexpresadas_2_tsv <-  paste(ruta_proteinas_transcriptas, "proteinas_coexpresadas_2.tsv", sep = "")
Proteinas_coexpresadas_3_1_tsv <-  paste(ruta_proteinas_transcriptas, "proteinas_coexpresadas_3_1.tsv", sep = "")
Proteinas_coexpresadas_3_2_tsv <-  paste(ruta_proteinas_transcriptas, "proteinas_coexpresadas_3_2.tsv", sep = "")

df_Proteinas_coexpresadas <- read_tsv(Proteinas_coexpresadas_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_Proteinas_coexpresadas_2 <- read_tsv(Proteinas_coexpresadas_2_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_Proteinas_coexpresadas_3_1 <- read_tsv(Proteinas_coexpresadas_3_1_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")
df_Proteinas_coexpresadas_3_2 <- read_tsv(Proteinas_coexpresadas_3_2_tsv, show_col_types = FALSE) %>%  filter(!is.na(product) & product != "")

#Filtrar proteínas en .fasta de las que son predichas comoe expresión del FTs

# Extraer los protein_id del archivo .tsv
proteins_ids <- list()

for(factor in factores_transcripcion){
  proteins_ids[[factor]] <- df_proteinas_tsv[[factor]]$protein_id
}

protein_coexpresadas_ids <- df_Proteinas_coexpresadas$protein_id
protein_coexpresadas_2_ids <- df_Proteinas_coexpresadas_2$protein_id
protein_coexpresadas_3_1_ids <- df_Proteinas_coexpresadas_3_1$protein_id
protein_coexpresadas_3_2_ids <- df_Proteinas_coexpresadas_3_2$protein_id

#  Leer el archivo FASTA con las proteínas
proteinas_fasta <- readAAStringSet("/home/Nephelim/Descargas/Beauveria bassiana-Genoma completo/ncbi_dataset/data/GCF_000280675.1/protein.faa")
nombres_proteinas_fasta <- names(proteinas_fasta)

# Filtrar las proteínas que están en el .tsv usando protein_id
# Convertir las cabeceras del FASTA (que contienen los protein_id) en un vector
fasta_ids <- sapply(strsplit(names(proteinas_fasta), " "), `[`, 1)

# Filtramos las proteínas cuyos protein_id están en el archivo .tsv
proteinas_filtradas <- list()

for(factor in factores_transcripcion){
  proteinas_filtradas[[factor]] <- proteinas_fasta[fasta_ids %in% proteins_ids[[factor]]]
}

proteinas_filtradas_coexpresadas <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_ids]
proteinas_filtradas_coexpresadas_2 <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_2_ids]
proteinas_filtradas_coexpresadas_3_1 <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_3_1_ids]
proteinas_filtradas_coexpresadas_3_2 <- proteinas_fasta[fasta_ids %in% protein_coexpresadas_3_2_ids]

# Guardar el resultado en un nuevo archivo FASTA
for(factor in factores_transcripcion){
  writeXStringSet(proteinas_filtradas[[factor]] , paste(ruta_proteinas_transcriptas, "proteinas_filtradas_", factor, ".fasta", sep = ""))
}

writeXStringSet(proteinas_filtradas_coexpresadas, paste(ruta_proteinas_transcriptas, "proteinas_filtradas_coexpresadas.fasta", sep =""))
writeXStringSet(proteinas_filtradas_coexpresadas_2, paste(ruta_proteinas_transcriptas, "proteinas_filtradas_coexpresadas_2.fasta", sep =""))
writeXStringSet(proteinas_filtradas_coexpresadas_3_1, paste(ruta_proteinas_transcriptas, "proteinas_filtradas_coexpresadas_3_1.fasta", sep =""))
writeXStringSet(proteinas_filtradas_coexpresadas_3_2, paste(ruta_proteinas_transcriptas, "proteinas_filtradas_coexpresadas_3_2.fasta", sep =""))

# Leer el archivo ignorando las líneas que comienzan con #
SignalP_data <- list()

ruta_signal_P <- ""

for(factor in factores_transcripcion){
  SignalP_data[[factor]] <- read_tsv(paste(ruta_signal_P, "Signal_P_6_", factor, ".tsv", sep = ""),
                                     comment = "#", show_col_types = FALSE, col_names = TRUE)
}

SignalP_data_ProteinasCoexpresadas <- read_tsv("/home/Nephelim/Descargas/Signal_P/Signal_P_ProteinasCoexpresadas.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_ProteinasCoexpresadas_2 <- read_tsv("/home/Nephelim/Descargas/Signal_P/Signal_P_ProteinasCoexpresadas_2.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_ProteinasCoexpresadas_3_1 <- read_tsv("/home/Nephelim/Descargas/Signal_P/Signal_P_ProteinasCoexpresadas_3_1.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)
SignalP_data_ProteinasCoexpresadas_3_2 <- read_tsv("/home/Nephelim/Descargas/Signal_P/Signal_P_ProteinasCoexpresadas_3_2.tsv",
                         comment = "#", show_col_types = FALSE, col_names = TRUE)

for(factor in factores_transcripcion){
  colnames(SignalP_data[[factor]]) <- c("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)","TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS_Position")
}

colnames(SignalP_data_ProteinasCoexpresadas) <- c("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)","TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS_Position")
colnames(SignalP_data_ProteinasCoexpresadas_2) <- c("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)","TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS_Position")
colnames(SignalP_data_ProteinasCoexpresadas_3_1) <- c("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)","TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS_Position")
colnames(SignalP_data_ProteinasCoexpresadas_3_2) <- c("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)","TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS_Position")

for(factor in factores_transcripcion){
  SignalP_data[[factor]] <- SignalP_data[[factor]] %>% filter(!is.na(ID) & ID != "")
  SignalP_data[[factor]] <- SignalP_data[[factor]] %>% mutate(ID = str_extract(ID, "^[^ ]+"))
  SignalP_data[[factor]] <- SignalP_data[[factor]] %>% filter(!is.na(CS_Position))
}

SignalP_data_ProteinasCoexpresadas <- SignalP_data_ProteinasCoexpresadas %>% filter(!is.na(ID) & ID != "")
SignalP_data_ProteinasCoexpresadas <- SignalP_data_ProteinasCoexpresadas %>% mutate(ID = str_extract(ID, "^[^ ]+"))
SignalP_data_ProteinasCoexpresadas <- SignalP_data_ProteinasCoexpresadas %>% filter(!is.na(CS_Position))

SignalP_data_ProteinasCoexpresadas_2 <- SignalP_data_ProteinasCoexpresadas_2 %>% filter(!is.na(ID) & ID != "")
SignalP_data_ProteinasCoexpresadas_2 <- SignalP_data_ProteinasCoexpresadas_2 %>% mutate(ID = str_extract(ID, "^[^ ]+"))
SignalP_data_ProteinasCoexpresadas_2 <- SignalP_data_ProteinasCoexpresadas_2 %>% filter(!is.na(CS_Position))

SignalP_data_ProteinasCoexpresadas_3_1 <- SignalP_data_ProteinasCoexpresadas_3_1 %>% filter(!is.na(ID) & ID != "")
SignalP_data_ProteinasCoexpresadas_3_1 <- SignalP_data_ProteinasCoexpresadas_3_1 %>% mutate(ID = str_extract(ID, "^[^ ]+"))
SignalP_data_ProteinasCoexpresadas_3_1 <- SignalP_data_ProteinasCoexpresadas_3_1 %>% filter(!is.na(CS_Position))

SignalP_data_ProteinasCoexpresadas_3_2 <- SignalP_data_ProteinasCoexpresadas_3_2 %>% filter(!is.na(ID) & ID != "")
SignalP_data_ProteinasCoexpresadas_3_2 <- SignalP_data_ProteinasCoexpresadas_3_2 %>% mutate(ID = str_extract(ID, "^[^ ]+"))
SignalP_data_ProteinasCoexpresadas_3_2 <- SignalP_data_ProteinasCoexpresadas_3_2 %>% filter(!is.na(CS_Position))

# Ver los primeros registros
head(SignalP_data[["BbAP1"]])
head(SignalP_data[["BbMsn2"]])
head(SignalP_data_ProteinasCoexpresadas)

# Pasar a data.frame
SignalP_data_df <- list()

for(factor in factores_transcripcion){
  SignalP_data_df[[factor]] <- as.data.frame(SignalP_data[[factor]])
}

SignalP_data_ProteinasCoexpresadas_df <- as.data.frame(SignalP_data_ProteinasCoexpresadas)
SignalP_data_ProteinasCoexpresadas_df_2 <- as.data.frame(SignalP_data_ProteinasCoexpresadas_2)
SignalP_data_ProteinasCoexpresadas_df_3_1 <- as.data.frame(SignalP_data_ProteinasCoexpresadas_3_1)
SignalP_data_ProteinasCoexpresadas_df_3_2 <- as.data.frame(SignalP_data_ProteinasCoexpresadas_3_2)

# Obtener los protein_id
proteinas_ids_SignalP <- list()

for(factor in factores_transcripcion){
  proteinas_ids_SignalP[[factor]] <- SignalP_data_df[[factor]]$ID
}

proteinas_id_ProteinasCoexpresadas_Signal_P <- SignalP_data_ProteinasCoexpresadas_df$ID
proteinas_id_ProteinasCoexpresadas_Signal_P_2 <- SignalP_data_ProteinasCoexpresadas_df_2$ID
proteinas_id_ProteinasCoexpresadas_Signal_P_3_1 <- SignalP_data_ProteinasCoexpresadas_df_3_1$ID
proteinas_id_ProteinasCoexpresadas_Signal_P_3_2 <- SignalP_data_ProteinasCoexpresadas_df_3_2$ID

# Filtrar las proteinas cuyo protein_id están en el archivo .tsv
proteinas_filtradas_SignalP <- list()

for(factor in factores_transcripcion){
  proteinas_filtradas_SignalP[[factor]] <- proteinas_fasta[fasta_ids %in% proteinas_ids_SignalP[[factor]]]
}

proteinas_filtradas_ProteinasCoexpresadas_SignalP <- proteinas_fasta[fasta_ids %in% proteinas_id_ProteinasCoexpresadas_Signal_P]
proteinas_filtradas_ProteinasCoexpresadas_SignalP_2 <- proteinas_fasta[fasta_ids %in% proteinas_id_ProteinasCoexpresadas_Signal_P_2]
proteinas_filtradas_ProteinasCoexpresadas_SignalP_3_1 <- proteinas_fasta[fasta_ids %in% proteinas_id_ProteinasCoexpresadas_Signal_P_3_1]
proteinas_filtradas_ProteinasCoexpresadas_SignalP_3_2 <- proteinas_fasta[fasta_ids %in% proteinas_id_ProteinasCoexpresadas_Signal_P_3_2]

#Guardar el resultado en un archivo .FASTA

for(factor in factores_transcripcion){
  writeXStringSet(proteinas_filtradas_SignalP[[factor]], paste(ruta_signal_P, "proteinas_filtradas_", factor, "_Signal_P.fasta", sep =""))
}

writeXStringSet(proteinas_filtradas_ProteinasCoexpresadas_SignalP, paste(ruta_signal_P, "proteinas_filtradas_", "ProteinasCoexpresadas", "_Signal_P.fasta", sep =""))
writeXStringSet(proteinas_filtradas_ProteinasCoexpresadas_SignalP_2, paste(ruta_signal_P, "proteinas_filtradas_", "ProteinasCoexpresadas_2", "_Signal_P.fasta", sep =""))
writeXStringSet(proteinas_filtradas_ProteinasCoexpresadas_SignalP_3_1, paste(ruta_signal_P, "proteinas_filtradas_", "ProteinasCoexpresadas_3_1", "_Signal_P.fasta", sep =""))
writeXStringSet(proteinas_filtradas_ProteinasCoexpresadas_SignalP_3_2, paste(ruta_signal_P, "proteinas_filtradas_", "ProteinasCoexpresadas_3_2", "_Signal_P.fasta", sep =""))

# Crear una lista vacía para almacenar los data frames filtrados
SignalP_annotation <- list()

for (factor in factores_transcripcion) {
  if (factor %in% names(anno_dfs) & factor %in% names(proteinas_ids_SignalP)) {
    SignalP_annotation[[factor]] <- anno_dfs[[factor]] %>% filter(protein_id %in% proteinas_ids_SignalP[[factor]])
  }
}

# Convertir la lista en un único data frame combinando los resultados de cada factor
SignalP_annotation_df <- bind_rows(SignalP_annotation, .id = "Factor")



#Contar el número de genes por factor en SignalP_annotation_df
df_genes_SignalP <- SignalP_annotation_df %>%
  group_by(Factor) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Factor = factor(Factor, levels = factores_transcripcion))  # Mantener el orden original

# Gráfico de barras
ggplot(df_genes_SignalP, aes(x = Factor, y = Count, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "plasma") +  # Colores llamativos
  labs(x = "Factor de Transcripción", y = "Número de genes en SignalP",
       title = "Comparación de Cantidad de Genes en SignalP por Factor de Transcripción") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 32, face = "bold", family = "Helvetica"),
        axis.text.x = element_text(family = "Helvetica", size = 25, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 25),
        axis.title.x = element_text(family = "Helvetica", size = 22, face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 22, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank()) +
  coord_flip()  # Orientación horizontal para mejor visualización


