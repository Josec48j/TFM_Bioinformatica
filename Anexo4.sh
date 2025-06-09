#!/bin/bash

# Definir rutas
INPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/3_FastQC-trimmed"
OUTPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/6_STAR_Aligned"
GENOME_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/5_STAR_Index"

# Número de hilos
THREADS=32

# Bucle para alinear todas las muestras
for R1 in "$INPUT_DIR"/*_R1_001_val_1.fq.gz; do
    base=$(basename "$R1" _R1_001_val_1.fq.gz)
    R2="$INPUT_DIR/${base}_R2_001_val_2.fq.gz"

    echo "Mapeando: $base"

    STAR --runThreadN "$THREADS" \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix "$OUTPUT_DIR/${base}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped None \
         --outSAMattributes NH HI NM MD AS

    echo "Finalizado: $base"
done

echo "🚀 Mapeo de todas las muestras completado!"

