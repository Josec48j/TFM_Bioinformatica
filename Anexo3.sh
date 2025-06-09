#!/bin/bash

# Definir rutas
GENOME_DIR="/home/alumno02/TFM/AnalisisDeRNAseq/Output/5_STAR_Index"
GENOME_FASTA="/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/GCF_900880675.1_fSpaAur1.1_genomic.fna"
GTF_FILE="/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/Annotations.gtf"

# Calcular sjdbOverhang: longitud de lectura - 1 (100 - 1 = 99)
SJDB_OVERHANG=99

# Ejecutar STAR en modo indexación con el ajuste recomendado
STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --genomeSAindexNbases 13
