# Herramientas-Tesis-Licenciatura-Biotecnologia
Herramientas de programación y scripts empleados en lenguaje R y Python para confeccionar el análisis de datos genómicos, proteícos y predicciones de expresión de proteinas y sus propiedades.

En la carpeta "Predicción proteinas por FTs" se encuentran los siguientes Scripts en R :

-"Conversion-pwm-a-meme.R" : Para la conversión de las matrices de interacción de un motivo proteico a ADN al formato .meme (de la Base de Datos JASPAR). Es un archivo de "position weight matrix" obtenido en la base de datos CIS-BP.

-"Analisis_tipo_ChIP_Seq_Simulacion_FIMO.R" : Análisis tipo ChIP-seq de los picos de interacción de los FTs al genoma de B. bassiana obtenidos por simulación del Servidor FIMO-MEME.

-"Algoritmos_filtrado_proteinas_potencialmente_transcritas.R" : Algoritmos y funciones para filtrar los genes que sean proteícos ordenandolos según su locus_tag en un data_frame y filtrar las proteínas que tenga alguna señal Signal_P obtenida por el servidor Signal P 6.0 DTU Health Tech
