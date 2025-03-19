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
  cat("Background letter frequencies\n", file=output_meme, append=TRUE)
  cat("A 0.25 C 0.25 G 0.25 T 0.25\n\n", file=output_meme, append=TRUE)
  cat(paste("MOTIF", motif_name, "\n\n"), file=output_meme, append=TRUE)
  cat(paste("letter-probability matrix: alength= 4 w=", w, "\n"), file=output_meme, append=TRUE)

  # Guardar la matriz PWM en formato MEME
  write.table(t(pwm_data), file=output_meme, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

# Uso de la función
convert_pwm_to_meme(nombreArchivopwm.tsv, nombreArchivo.meme, motif_name="nombre del motivo")
