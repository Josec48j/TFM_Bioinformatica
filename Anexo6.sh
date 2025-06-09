#!/bin/bash

# Directorios
GTF_FILE="/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/Annotations.gtf"
INPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/6_STAR_Aligned/"
OUTPUT_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/7_FeatureCounts"

# Número de hilos
THREADS=32

# Ejecutar featureCounts en todos los archivos BAM
featureCounts -T $THREADS \
  -a "$GTF_FILE" \
  -o "$OUTPUT_DIR/gene_counts.txt" \
  -s 2 \
  -p --countReadPairs \
  -g gene_id \
  -t gene \
  "$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam

echo "🚀 Conteo de lecturas completado!"

