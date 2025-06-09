#!/bin/bash

# Directorios y archivo de anotaciones
INPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/6_STAR_Aligned"
OUTPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/6_STAR_Aligned/Qualimap"
GTF_FILE="/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/Annotations.gtf"

# Crear el directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Recorrer todos los archivos BAM en el directorio de entrada
for bam in "$INPUT_DIR"/*.bam; do
    base=$(basename "$bam" .bam)
    echo "Evaluando la calidad para $base..."
    
    qualimap rnaseq \
      -bam "$bam" \
      -gtf "$GTF_FILE" \
      -outdir "$OUTPUT_DIR/$base" \
      -outfile "${base}_report.pdf" \
      -outformat PDF \
      --java-mem-size=24G

    echo "Calidad evaluada para $base"
done

echo "🚀 Análisis de calidad completado!"
