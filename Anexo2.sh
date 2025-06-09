#!/bin/bash

# Directorios
INPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/data/1_Archivos_fastq/"
OUTPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/3_FastQC-trimmed"

# Bucle para procesar todos los archivos R1 y R2 en el directorio de entrada
for file in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    base=$(basename "$file" _R1_001.fastq.gz)
    echo "Procesando: $base"

    R1="$INPUT_DIR/${base}_R1_001.fastq.gz"
    R2="$INPUT_DIR/${base}_R2_001.fastq.gz"

 #  Trim Galore (para quitar adaptadores y correr FastQC en modo paired)
    trim_galore --paired --fastqc -o "$OUTPUT_DIR" "$R1" "$R2"
    echo "Finalizado: $base"
done

echo "🚀 Trimming de todas las muestras paired-end completado!"
