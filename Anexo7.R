# ==== Librerías para manipulación de datos ====
library(tibble)       # Tibbles y funciones auxiliares para marcos de datos
library(dplyr)        # Gramática de datos: filter, mutate, select, summarize, join…
library(tidyr)        # Transformaciones wide ↔ long: pivot_longer, pivot_wider, etc.
library(purrr)        # Programación funcional: map, reduce, walk sobre listas y vectores
library(stringr)      # Operaciones con texto: detección, reemplazo y cambio de caso
library(forcats)      # Manejo de factores: reordenar, recodificar niveles, etc.

# ==== Librerías para I/O de Excel ====
library(openxlsx)     # Leer y escribir archivos Excel con múltiples hojas
library(readxl)       # Leer archivos Excel (xls, xlsx)

# ==== Librerías para control de calidad y normalización ====
library(NOISeq)       # QC de RNA-seq: bias plots, addData, explo.plot, etc.
library(edgeR)        # Conteo por millón (cpm) y normalización TMM
library(sva)          # Corrección de batch-effect con ComBat

# ==== Librerías para análisis de expresión diferencial ====
library(DESeq2)       # Modelo DESeq2 para análisis de expresión diferencial
library(apeglm)       # Shrinkage de log2FC para DESeq2 (lfcShrink)

# ==== Librerías para anotaciones genómicas ====
library(rtracklayer)  # Importar GTF/BED para extraer biotipos y coordenadas
library(AnnotationDbi)# Mapear entre IDs (ENSEMBL ↔ ENTREZ, etc.)
library(org.Dr.eg.db) # Base de datos de anotación de Danio rerio

# ==== Librerías para visualización ====
library(ggplot2)      # Gráficos basados en la gramática de gráficos
library(ggrepel)      # Etiquetas no solapadas en gráficos

# ==== Librerías para visualización de DEGs ====
library(EnhancedVolcano)  # Volcanes enriquecidos para resultados de DESeq2

# ==== Librerías para visualización de rutas metabólicas ====
library(pathview)         # Visualización de rutas KEGG con datos genómicos

# ==== Librerías para enriquecimiento funcional ====
library(clusterProfiler)  # Análisis ORA y GSEA sobre GO y KEGG
library(enrichplot)       # Gráficos para resultados de enriquecimiento
library(ReactomePA)       # Enriquecimiento de rutas específicas de Reactome

# ==== Librerías para heatmaps avanzados ====
library(ComplexHeatmap)   # Heatmaps con anotaciones y rasterizado
library(circlize)         # Paletas de color y colorRamp2 para heatmaps


###############################################################################
# I. Definir rutas y cargar datos
###############################################################################
# 1. rutas comunes
annot_file <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/annotation_dorada.tsv"
outdir     <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/1_Control_calidad"
outdir2    <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/2_DEG"
outdir3    <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/3_Enriquecimiento"
orth_file <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/3_Ortologos_Zebra_con_Entrez_GO.xlsx"
Ortologos_file   <- read_excel("/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/2_Ortologos_pez_zebra_con_GO.xlsx")
gtf        <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/3_Genome_and_annotations/Annotations.gtf"

# 2. Lista de tejidos y sus rutas/“letras”
tissues <- list(
  Piel     = list(count_file = "/home/alumno02/TFM/AnalisisDeRNAseq/Output/7_FeatureCounts/Piel/gene_counts.txt",
                  letter     = "P"),
  Cerebro  = list(count_file = "/home/alumno02/TFM/AnalisisDeRNAseq/Output/7_FeatureCounts/Cerebro/gene_counts.txt",
                  letter     = "B")
)

# 3. Función genérica de carga + metadatos
process_counts <- function(count_file, letter) {
  # Leer y filtrar
  fc     <- read.table(count_file, header = TRUE, row.names = 1,
                       sep = "\t", stringsAsFactors = FALSE)
  counts <- fc[, -(1:5)]
  counts <- counts[rowSums(counts) > 0, ]
  rownames(counts) <- sub("\\..*$", "", rownames(counts))
  
  # Extraer sample_ids
  pattern     <- paste0("^.*?(\\d+", letter, ")_S\\d+.*$")
  sample_ids  <- gsub(pattern, "\\1", colnames(counts))
  
  # Función de agrupamiento
  get_group <- function(s) {
    if      (s %in% paste0(1:4,  letter)) "Control"
    else if (s %in% paste0(7:10, letter)) "Carragenina"
    else if (s %in% paste0(13:16,letter)) "SST6"
    else if (s %in% paste0(25:28,letter)) "Carragenina_SST6"
    else NA
  }
  
  # Crear coldata
  groups  <- factor(sapply(sample_ids, get_group),
                    levels = c("Control","Carragenina","SST6","Carragenina_SST6"))
  coldata <- data.frame(condition = groups, row.names = sample_ids)
  colnames(counts) <- sample_ids
  
  # Devolver resultados
  list(counts = counts, coldata = coldata)
}

# 4. Ejecutar para cada tejido y asignar a tu espacio de trabajo

results <- lapply(tissues, function(x) {
  process_counts(count_file = x$count_file, letter = x$letter)
})

# 5. para traerlos al workspace global
for (t in names(results)) {
  assign(paste0("counts_",    t), results[[t]]$counts)
  assign(paste0("coldata_",   t), results[[t]]$coldata)
}

###############################################################################
# II. Control de calidad con NOISeq para múltiples tejidos
###############################################################################
# Leer anotaciones
annot <- read.delim(annot_file, stringsAsFactors = FALSE)  # gene_id | length | gc
annot$gene_id   <- sub("\\..*$", "", annot$gene_id)
rownames(annot) <- annot$gene_id

# Definir función genérica de control de calidad
qc_noiseq <- function(counts, coldata, tissue, outdir) {
  # Filtrar genes comunes entre counts y annot
  common <- intersect(rownames(counts), rownames(annot))
  counts_f <- counts[common, ]
  annot_f  <- annot[common, ]
  
  # Crear objeto NOISeq
  myData <- readData(data    = counts_f,
                     factors = coldata)
# 1) SESGOS DE LONGITUD Y GC
  # 2.3) Añadir longitud y GC
  length_vec <- setNames(annot_f$length, common)
  gc_vec     <- setNames(annot_f$gc,     common)
  myData <- addData(myData, length = length_vec, gc = gc_vec)
  
  # Plots de sesgos
  # Longitud
  lb_raw <- dat(myData, factor = "condition", type = "lengthbias", norm = FALSE)
  png(file.path(outdir, paste0("NOISeq_bias_length_raw_", tissue, ".png")),
      width = 2000, height = 2000, res = 300)
  explo.plot(lb_raw, samples = NULL, toplot = "global")
  dev.off()
  
  # GC
  gc_raw <- dat(myData, factor = "condition", type = "GCbias", norm = FALSE)
  png(file.path(outdir, paste0("NOISeq_bias_GC_raw_", tissue, ".png")),
      width = 2000, height = 2000, res = 300)
  explo.plot(gc_raw, samples = NULL, toplot = "global")
  dev.off()
  
  message("✔ QC NOISeq completado para ", tissue)
}

# Ejecutar para Piel y Cerebro
qc_noiseq(counts_Piel,    coldata_Piel,    tissue = "Piel",    outdir)
qc_noiseq(counts_Cerebro, coldata_Cerebro, tissue = "Cerebro", outdir)

# Lista con los datos de cada tejido
tissues <- list(
  Piel    = list(counts = counts_Piel,    coldata = coldata_Piel,    letter = "P"),
  Cerebro = list(counts = counts_Cerebro, coldata = coldata_Cerebro, letter = "B")
)

# 2) FILTRADO DE BAJO RECUENTO (CPM >1 en ≥2 muestras)
counts_norm <- list()
for(t in names(tissues)) {
  ct  <- tissues[[t]]$counts
  cpm_mat <- cpm(ct)
  keep    <- rowSums(cpm_mat > 1) >= 2
  message(t, ": retenidos ", sum(keep), " genes tras low-count filter")
  counts_norm[[t]] <- round(ct[keep, ])
}

# 3) CORRECCIÓN DE EFECTO DE LOTE (ComBat) — sólo si hay columna “batch”
counts_corr <- counts_norm
for(t in names(tissues)) {
  cd <- tissues[[t]]$coldata
  if("batch" %in% colnames(cd)) {
    mod  <- model.matrix(~ condition, cd)
    mat  <- as.matrix(counts_norm[[t]])
    cb   <- ComBat(dat = mat, batch = cd$batch, mod = mod)
    counts_corr[[t]] <- round(cb)
    message(t, ": batch corrected")
  } else {
    message(t, ": no hay columna 'batch', se salta ComBat")
  }
}

# 4) PCA exploratorio
for(t in names(tissues)) {
  cd  <- tissues[[t]]$coldata
  cn  <- counts_corr[[t]]
  
  log2_mat <- log2(cn + 1)
  pca      <- prcomp(t(log2_mat), center=TRUE, scale.=FALSE)
  pct1     <- round(100 * pca$sdev[1]^2 / sum(pca$sdev^2))
  pct2     <- round(100 * pca$sdev[2]^2 / sum(pca$sdev^2))
  
  pca_df <- data.frame(
    PC1       = pca$x[,1],
    PC2       = pca$x[,2],
    sample    = rownames(pca$x),
    condition = cd[rownames(pca$x),"condition"],
    stringsAsFactors = FALSE
  )
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, fill = condition)) +
    geom_point(shape = 21, colour = "black", size = 5, stroke = 1) +
    geom_text(aes(label = sample), vjust = -1.5, fontface = "bold") +
    scale_fill_manual(values = c(
      Control           = "#8dd3c7",
      Carragenina       = "#fb8072",
      SST6              = "#80b1d3",
      Carragenina_SST6  = "#bebada"
    )) +
    labs(
      title = paste0(t, " – PCA"),
      x     = paste0("PC1 (", pct1, "%)"),
      y     = paste0("PC2 (", pct2, "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.title       = element_blank(),
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.box.margin  = margin(0, 0, 0, 12),
      plot.title         = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0(t, "_PCA.png")),
    plot   = p_pca,
    width  = 7,
    height = 7,
    dpi    = 300
  )
}

# 5) MAD OUTLIER DETECTION — identificar muestras atípicas
madScore <- function(mat, sheet, groupby="condition", thr=-5) {
  mat <- mat[, rownames(sheet)]
  mats <- lapply(split(rownames(sheet), sheet[[groupby]]),
                 function(s) mat[, s, drop=FALSE])
  mads <- do.call(rbind, lapply(names(mats), function(g) {
    M <- mats[[g]]; C <- cor(M)
    corr_in   <- sapply(seq_len(ncol(C)), function(i) mean(C[i,-i]))
    corr_noin <- mean(C[upper.tri(C)])
    diffs     <- corr_in - corr_noin
    md        <- median(diffs)
    madv      <- median(abs(diffs - md))
    data.frame(group = g,
               mad   = (diffs - md)/(madv * 1.4826),
               row.names = colnames(M),
               check.names = FALSE)
  }))
  mads$outlier <- mads$mad < thr
  mads
}

for(t in names(tissues)) {
  cd   <- tissues[[t]]$coldata
  mat  <- log2(counts_corr[[t]] + 1)
  pd   <- madScore(mat, cd, "condition", -5)
  outliers <- rownames(pd)[pd$outlier]
  
  # Plot con fondo blanco, ejes y marco
  df <- pd %>% rownames_to_column("sample")
  p_mad <- ggplot(df, aes(group, mad, label=sample)) +
    geom_point(aes(color = outlier), size = 3) +
    geom_text(hjust = -0.2, size = 3, alpha = 0.6) +
    geom_hline(yintercept = -5, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "orange")) +
    labs(title = paste0(t, " – MAD outliers"),
         x     = "Grupo",
         y     = "MAD score",
         color = "Outlier") +
    theme_classic(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      axis.title    = element_text(face = "bold"),
      axis.text     = element_text(face = "bold"),
      panel.border  = element_rect(colour = "black", fill = NA),
      legend.position = "right"
    )
  
  ggsave(
    file.path(outdir, paste0(t, "_MAD_outliers.png")),
    plot  = p_mad,
    width = 6, height = 4, dpi = 300
  )
  
  # Eliminar outliers de counts y coldata
  if (length(outliers) > 0) {
    counts_corr[[t]] <- counts_corr[[t]][, !colnames(counts_corr[[t]]) %in% outliers]
    tissues[[t]]$coldata <- tissues[[t]]$coldata[!rownames(tissues[[t]]$coldata) %in% outliers, , drop = FALSE]
    message(t, ": eliminadas muestras outlier:", paste(outliers, collapse = ", "))
  } else {
    message(t, ": no se detectaron outliers")
  }
}

# Creamos copias de los contajes sin outliers:
# Para Piel
counts_Piel_SinOutliers  <- counts_corr[["Piel"]]
coldata_Piel_SinOutliers <- tissues[["Piel"]]$coldata

# Para Cerebro
counts_Cerebro_SinOutliers  <- counts_corr[["Cerebro"]]
coldata_Cerebro_SinOutliers <- tissues[["Cerebro"]]$coldata

###############################################################################
# III. Expresión diferencial
###############################################################################
# 1) Preparamos la lista de datos para DE, incluyendo los datos sin la muestra 3
DE_samples <- list(
  Piel   = list(counts  = counts_Piel_SinOutliers,
                coldata = coldata_Piel_SinOutliers),
  Cerebro= list(counts  = counts_Cerebro_SinOutliers,
                coldata = coldata_Cerebro_SinOutliers)
)

# 2) Construcción y ajuste del modelo DESeq2 para ambos conjuntos
for(t in names(DE_samples)) {
  dat <- DE_samples[[t]]
  
  # Creamos el DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = dat$counts,
    colData   = dat$coldata,
    design    = ~ condition
  )
  # Filtramos genes de sum(counts) >= 1
  dds <- dds[rowSums(counts(dds)) >= 1, ]
  
  # Ajustamos DESeq
  dds <- DESeq(dds, parallel = TRUE)
  
  # Guardamos el objeto en el workspace
  assign(paste0("dds_", t), dds)
  message("✔ DESeq2 ajustado para ", t)
}

# 3) PCA SOBRE CUENTAS NORMALIZADAS POR DESeq2 para ambos conjuntos
for(t in names(DE_samples)) {
  # Obtener el objeto DESeq ajustado
  dds   <- get(paste0("dds_", t))
  # Metadatos para el conjunto actual
  cd    <- DE_samples[[t]]$coldata
  
  # Extraer matriz normalizada y hacer PCA
  mat_log <- log2(counts(dds, normalized = TRUE) + 1)
  pca     <- prcomp(t(mat_log), center = TRUE, scale. = FALSE)
  pct     <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  # Preparar dataframe para ggplot
  pca_df <- data.frame(
    PC1       = pca$x[,1],
    PC2       = pca$x[,2],
    sample    = rownames(pca$x),
    condition = cd[rownames(pca$x), "condition"],
    stringsAsFactors = FALSE
  )
  
  # Crear el gráfico
  p_pca <- ggplot(pca_df, aes(PC1, PC2, fill = condition)) +
    geom_point(shape = 21, colour = "black", size = 5, stroke = 1) +
    geom_text(aes(label = sample), vjust = -1.5, fontface = "bold") +
    scale_fill_manual(values = c(
      Control           = "#8dd3c7",
      Carragenina       = "#fb8072",
      SST6              = "#80b1d3",
      Carragenina_SST6  = "#bebada"
    )) +
    labs(
      title = paste0(t, " – PCA DESeq2"),
      x     = paste0("PC1 (", pct[1], "%)"),
      y     = paste0("PC2 (", pct[2], "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major  = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(face = "bold"),
      legend.title      = element_blank(),
      legend.position   = "right",
      legend.background = element_blank(),
      plot.title        = element_text(face = "bold", hjust = 0.5)
    )
  
  # Guardar en outdir2
  ggsave(
    filename = file.path(outdir2, paste0(t, "_PCA_DESeq2.png")),
    plot     = p_pca,
    width    = 7,
    height   = 7,
    dpi      = 300
  )
}

# 4) EXTRAER RESULTADOS y GUARDAR EN EXCEL (hoja por contraste)
# Definir pares de comparación a partir de los datos “sin outliers”
conds <- levels(DE_samples[["Piel"]]$coldata$condition)
pairs <- combn(conds, 2, simplify = FALSE)
names(pairs) <- sapply(pairs, function(x) paste0(x[1], "_vs_", x[2]))

# Iterar sobre cada tejido
for (t in names(DE_samples)) {
  # Trabajamos con el objeto dds ajustado para este tejido
  dds <- get(paste0("dds_", t))
  
  # Creamos un nuevo workbook
  wb <- createWorkbook()
  
  # Recorremos cada contraste y volcamos los resultados en una hoja
  for (nm in names(pairs)) {
    a <- pairs[[nm]][1]
    b <- pairs[[nm]][2]
    
    res_df <- results(dds, contrast = c("condition", a, b)) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      filter(!is.na(padj)) %>%
      arrange(desc(log2FoldChange))
    
    addWorksheet(wb, nm)
    writeData(wb, nm, res_df)
  }
  
  # Guardar el Excel directamente en outdir2, con nombre por tejido
  saveWorkbook(
    wb,
    file      = file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx")),
    overwrite = TRUE
  )
  message("✔ Excel DESeq2 guardado para ", t,
          " en ", file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx")))
}

# 5) Resumen rápido de DEGs por contraste
for (t in names(DE_samples)) {
  
  excel_file <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  
  deg_summary <- lapply(names(pairs), function(nm){
    df <- read.xlsx(excel_file, sheet = nm)
    tibble(
      Contrast = nm,
      Up   = sum(df$padj < 0.05 & df$log2FoldChange >  1),
      Down = sum(df$padj < 0.05 & df$log2FoldChange < -1)
    )
  }) |>
    bind_rows()
  
  write.xlsx(
    deg_summary,
    file = file.path(outdir2, paste0(t, "_DEG_summary.xlsx")),
    overwrite = TRUE
  )
  message("✔ Tabla de resumen guardada para ", t)
}

# 6) VOLCANO PLOTS (EnhancedVolcano)
# Pares ya definidos arriba
for (t in names(DE_samples)) {
  # Ruta al Excel de resultados DE para este tejido
  excel_file <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  
  for (nm in names(pairs)) {
    # Leer resultados del contraste
    res_df <- read.xlsx(excel_file, sheet = nm)
    
    # Generar volcano plot
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab            = res_df$gene_id,
      x              = "log2FoldChange",
      y              = "padj",
      title          = paste0(t, " | ", nm),
      pCutoff        = 0.05,
      FCcutoff       = 1.0,
      pointSize      = 3.0,
      labSize        = 3.0,
      legendPosition = "right",
      legendLabSize  = 10,
      legendIconSize = 3.0,
      col            = c("gray70","forestgreen","royalblue","red2"),
      cutoffLineCol  = "black",
      cutoffLineType = "dashed",
      caption        = "Fuente: DESeq2 + EnhancedVolcano",
      xlab           = bquote(~Log[2]~FC),
      ylab           = bquote(-Log[10]~adj~p)
    )
    
    # Guardar PNG
    ggsave(
      filename = file.path(outdir2, paste0(t, "_EnhancedVolcano_", nm, ".png")),
      plot     = volcano_plot,
      width    = 8,
      height   = 6,
      dpi      = 300
    )
  }
  message("✔ Volcano plots creados para ", t)
}

# 7) Heat-map por condición
for (t in names(DE_samples)) {
  # VST
  dds <- get(paste0("dds_", t))
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  # Listado de DEGs
  deg_list <- lapply(names(pairs), function(nm) {
    read.xlsx(
      file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx")),
      sheet = nm
    ) %>%
      filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
      pull(gene_id)
  })
  deg_genes <- unique(unlist(deg_list))
  if (length(deg_genes) < 2) next
  
  # Matriz de medias por condición
  mat <- assay(vsd)[deg_genes, , drop = FALSE]
  conds <- colData(dds)$condition
  group_means <- sapply(levels(conds), function(cn)
    rowMeans(mat[, conds == cn, drop = FALSE])
  )
  colnames(group_means) <- levels(conds)
  
  # Z-score por gen
  mat_scl <- t(scale(t(group_means)))
  
  # Heatmap
  ht <- Heatmap(
    mat_scl,
    name             = "Z-score",
    show_row_names   = FALSE,
    cluster_rows     = TRUE,
    cluster_columns  = TRUE,
    column_names_rot = 45,
    column_names_gp  = gpar(fontsize = 10),
    column_title     = paste0("Heatmap DEGs por grupo - ", t),
    use_raster       = TRUE
  )
  
  # Guardar PNG
  png(
    filename = file.path(outdir2, paste0(t, "_Heatmap_DEGs_groups.png")),
    width    = 6, height = 5, units = "in",
    res      = 300, type = "cairo"
  )
  draw(ht)
  dev.off()
}

# 7) Boxplot de switchers (genes que invierten dirección con SST-6) + Excel
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)

switcher_boxplot <- function(tissue, outdir = outdir2){
  xls <- file.path(outdir, paste0(tissue,"_DESeq2_all_comparisons.xlsx"))
  
  ## Genes switcher
  carr   <- read.xlsx(xls, sheet = "Control_vs_Carragenina")
  carrS6 <- read.xlsx(xls, sheet = "Control_vs_Carragenina_SST6")
  sw <- merge(carr, carrS6, by = "gene_id",
              suffixes = c(".Carr", ".CarrS6")) |>
    filter(padj.Carr  < .05,
           padj.CarrS6 < .05,
           sign(log2FoldChange.Carr) != sign(log2FoldChange.CarrS6))
  genes_sw <- sw$gene_id
  if (!length(genes_sw)){
    message(tissue, ": no hay switchers"); return(invisible())
  }
  
  ## Contrastes
  keep <- c("Control_vs_Carragenina",
            "Control_vs_SST6",
            "Control_vs_Carragenina_SST6",
            "Carragenina_vs_SST6",
            "Carragenina_vs_Carragenina_SST6")
  lfc_list <- lapply(keep, function(sh){
    df <- read.xlsx(xls, sheet = sh)
    sub <- df[df$gene_id %in% genes_sw,
              c("gene_id","log2FoldChange")]
    names(sub)[2] <- sh
    sub
  })
  m_df <- Reduce(function(x,y) merge(x,y,by="gene_id",all=TRUE),
                 lfc_list)
  
  ## Formato long
  long <- pivot_longer(m_df, -gene_id,
                       names_to="contrast", values_to="log2FC") |>
    mutate(contrast = factor(contrast, levels = keep))
  
  ## colores
  col_map <- c(
    Control_vs_Carragenina           = "#fb8072",
    Control_vs_SST6                  = "#80b1d3",
    Control_vs_Carragenina_SST6      = "#bebada",
    Carragenina_vs_SST6              = "#BE98A3",
    Carragenina_vs_Carragenina_SST6  = "#DD9DA6"
  )
  # Bordes más oscuros
  border_map <- col_map |> 
    lapply(function(cl){
      rgb(t(col2rgb(cl))/255 * 0.7)
    }) |> 
    unlist()
  
  ## Boxplot + caps + puntos
  p <- ggplot(long, aes(contrast, log2FC,
                        fill   = contrast,
                        colour = contrast)) +
    # caja sin outliers
    geom_boxplot(
      width         = 0.6,
      size          = 0.5,
      outlier.shape = NA
    ) +
    # Stat_boxplot se encargará de dibujar whiskers + caps
    stat_boxplot(
      geom  = "errorbar",
      width = 0.15,
      size  = 0.5
    ) +
    # Puntos individuales
    geom_jitter(
      width  = 0.15,
      size   = 1.3,
      alpha  = 0.8,
      colour = "black"
    ) +
    scale_fill_manual(values  = col_map) +
    scale_colour_manual(values = border_map) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      panel.grid      = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = paste0(tissue, " – Log2FC de genes switcher"),
      x     = NULL,
      y     = "log2(FC)"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(tissue,"_switchers_boxplot.png")),
    plot     = p,
    width    = 7, height = 5, dpi = 300
  )
  
  ## Excel con dos pestañas
  wb <- createWorkbook()
  addWorksheet(wb, "Switchers_only")
  writeData(wb, "Switchers_only",
            data.frame(gene_id = genes_sw),
            colNames = TRUE)
  addWorksheet(wb, "Log2FC_matrix")
  writeData(wb, "Log2FC_matrix", m_df, colNames = TRUE)
  saveWorkbook(wb,
               file = file.path(outdir,
                                paste0(tissue, "_Switchers_log2FC.xlsx")),
               overwrite = TRUE)
  
  message(tissue, ": boxplot + Excel con switchers generados")
}

# Ejemplo de uso:
switcher_boxplot("Piel")
switcher_boxplot("Cerebro")


# 7) Diagramas de Venn de DEGs Control vs tratamientos  +  Excel
library(VennDiagram)
library(openxlsx)
library(grid)

make_venn_vd <- function(tissue, outdir = outdir2) {
  # Ruta de tu Excel
  xls <- file.path(outdir2, paste0(tissue, "_DESeq2_all_comparisons.xlsx"))
  
  # Hoja → nombre de conjunto
  sheet_map <- c(
    Carragenina        = "Control_vs_Carragenina",
    SST6               = "Control_vs_SST6",
    `Carragenina+SST6` = "Control_vs_Carragenina_SST6"
  )
  
  # Extraemos listas de genes
  gene_sets <- lapply(sheet_map, function(sh) {
    df <- read.xlsx(xls, sheet = sh)
    subset(df, padj < 0.05 & abs(log2FoldChange) > 1)$gene_id
  })
  names(gene_sets) <- names(sheet_map)
  
  # Conteos de cada región
  a1   <- length(gene_sets$Carragenina)
  a2   <- length(gene_sets$SST6)
  a3   <- length(gene_sets$`Carragenina+SST6`)
  n12  <- length(intersect(gene_sets$Carragenina, gene_sets$SST6))
  n23  <- length(intersect(gene_sets$SST6, gene_sets$`Carragenina+SST6`))
  n13  <- length(intersect(gene_sets$Carragenina, gene_sets$`Carragenina+SST6`))
  n123 <- length(Reduce(intersect, gene_sets))
  
  # Colores base (sin alfa en el hex)
  # fill_cols <- c("#440154", "#21908D", "#FDE725")
  fill_cols <- c("#fb8072", "#80b1d3", "#825F97")

  # Dibujamos con bordes del mismo color, 40% de opacidad y sin decimales
  venn.plot <- draw.triple.venn(
    area1      = a1,
    area2      = a2,
    area3      = a3,
    n12        = n12,
    n23        = n23,
    n13        = n13,
    n123       = n123,
    category   = names(sheet_map),
    fill       = fill_cols,
    alpha      = rep(0.4, 3),   # 40% transparencia
    col        = fill_cols,     # bordes mismo color
    lwd        = 2,             # grosor de borde
    cex        = 1.3,           # tamaño de números
    cat.cex    = 1.2,           # tamaño de nombres
    cat.col    = fill_cols,     # color de nombres
    print.mode = c("raw", "percent"),
    sigdigs    = 0,             # sin decimales en los porcentajes
    euler.d    = FALSE,
    scaled     = FALSE
  )
  
  # Guardamos el gráfico
  png(
    filename = file.path(outdir, paste0(tissue, "_Venn_DEGs_vs_Control.png")),
    width    = 7, height = 6, units = "in", res = 300
  )
  grid.draw(venn.plot)
  dev.off()
  
  # Excel con listas
  carr  <- gene_sets$Carragenina
  sst6  <- gene_sets$SST6
  carrS <- gene_sets$`Carragenina+SST6`
  reg <- list(
    Exclusive_Carragenina      = setdiff(carr,  union(sst6, carrS)),
    Exclusive_SST6             = setdiff(sst6, union(carr,  carrS)),
    Exclusive_Carragenina_SST6 = setdiff(carrS, union(carr,  sst6)),
    Carragenina_and_SST6_only  = setdiff(intersect(carr, sst6), carrS),
    Carragenina_and_CarrS_only = setdiff(intersect(carr, carrS), sst6),
    SST6_and_CarrS_only        = setdiff(intersect(sst6, carrS), carr),
    Triple_Intersection        = Reduce(intersect, gene_sets)
  )
  wb <- createWorkbook()
  for (nm in names(reg)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, data.frame(Gene_ID = reg[[nm]]), colNames = FALSE)
  }
  saveWorkbook(
    wb,
    file.path(outdir, paste0(tissue, "_Venn_gene_lists.xlsx")),
    overwrite = TRUE
  )
}

# Ejecuta para todos los samples
invisible(lapply(names(DE_samples), make_venn_vd))


###############################################################################
# IV. ENRIQUECIMIENTO GO – PARA PIEL y CEREBRO
###############################################################################
# Mapeo Dorada_Symbol → Zebra_Entrez
df_orth <- readxl::read_excel(orth_file)
mapping <- unique(df_orth[, c("Dorada_Symbol","Zebra_Entrez")])

# Contrastres
conds <- levels(DE_samples[["Piel"]]$coldata$condition)
pairs <- combn(conds,2,simplify=FALSE)
names(pairs) <- sapply(pairs, function(x) paste0(x[1], "_vs_", x[2]))

# Carpetas por ontología
ont_dirs <- c(
  BP = "1_Biological_process",
  CC = "2_Cellular_Compartment",
  MF = "3_Molecular_Function"
)
for(d in ont_dirs) dir.create(file.path(outdir3,d), recursive=TRUE, showWarnings=FALSE)

# Loop por tejido y contraste
for (t in names(DE_samples)) {
  de_file <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  
  for (nm in names(pairs)) {
    # Creamos el workbook para este contraste
    wb_go <- openxlsx::createWorkbook()
    
    res_df <- openxlsx::read.xlsx(de_file, sheet = nm)
    names(res_df)[1] <- "Dorada_Symbol"
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # join a ENTREZ
    res_ent <- merge(res_df, mapping, by = "Dorada_Symbol")
    res_ent <- res_ent[!is.na(res_ent$Zebra_Entrez), ]
    if (nrow(res_ent) == 0) next
    
    # prerank
    geneList <- setNames(
      res_ent$log2FoldChange,
      as.character(res_ent$Zebra_Entrez)
    )
    geneList <- sort(geneList, decreasing = TRUE)
    
    for (ont in names(ont_dirs)) {
      ego <- clusterProfiler::gseGO(
        geneList     = geneList,
        OrgDb        = org.Dr.eg.db,
        keyType      = "ENTREZID",
        ont          = ont,
        minGSSize    = 10,
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )
      ego_df <- as.data.frame(ego) %>%
        mutate(
          Direction = ifelse(NES >= 0, "Activated", "Repressed")
        )
      if (nrow(ego_df) == 0) next
      
      # Count y GeneRatio
      ego_df$Count     <- lengths(strsplit(ego_df$core_enrichment, "/"))
      ego_df$GeneRatio <- ego_df$Count / ego_df$setSize
      
      # Top-10 y etiquetado (fijo los niveles de Direction)
      topn <- ego_df %>%
        arrange(p.adjust) %>%
        slice_head(n = 10) %>%
        mutate(
          Description = str_to_sentence(str_replace_all(Description, "-", " ")),
          Direction   = factor(
            ifelse(NES >= 0, "Activated", "Repressed"),
            levels = c("Activated","Repressed")
          )
        )
      
      # A) Si faltan Activados, añadimos un dummy row
      if (! "Activated" %in% topn$Direction) {
        topn <- dplyr::bind_rows(
          topn,
          data.frame(
            Description = NA,
            GeneRatio   = NA,
            Count       = NA,
            Direction   = factor("Activated", levels = c("Activated","Repressed"))
          )
        )
      }
      # B) Si faltan Reprimidos, idem
      if (! "Repressed" %in% topn$Direction) {
        topn <- dplyr::bind_rows(
          topn,
          data.frame(
            Description = NA,
            GeneRatio   = NA,
            Count       = NA,
            Direction   = factor("Repressed", levels = c("Activated","Repressed"))
          )
        )
      }
      
      # Bubble plot
      p <- ggplot(topn, aes(
        x     = GeneRatio,
        y     = factor(Description, levels = rev(Description)),
        size  = Count,
        fill  = Direction
      )) +
        geom_point(shape = 21, colour = "black", stroke = 0.5) +
        scale_fill_manual(
          name   = "Direction",
          values = c(Activated = "#fb8072", Repressed = "#80b1d3"),
          breaks = c("Activated", "Repressed"),
          drop   = FALSE
        ) +
        guides(fill = guide_legend(
          override.aes = list(
            size   = 5,
            shape  = 21,
            colour = "black",
            fill   = c("#fb8072", "#80b1d3")
          ),
          order          = 1,
          title.position = "top"
        )) +
        scale_size_continuous(range = c(3, 8)) +
        facet_grid(
          . ~ Direction,
          scales = "free_x",
          drop   = FALSE
        ) +
        scale_x_continuous(
          expand = expansion(mult = c(0.06, 0.06))
        ) +
        scale_y_discrete(na.translate = FALSE) +
        coord_cartesian(clip = "off") +
        labs(
          title = paste(t, nm, ont, sep = " | "),
          x     = "Gene ratio",
          y     = NULL
        ) +
        theme_bw(base_size = 14) +
        theme(
          panel.spacing       = unit(1, "lines"),
          panel.background    = element_rect(fill = "white"),
          panel.border        = element_rect(color = "black"),
          panel.grid.major    = element_blank(),
          panel.grid.minor    = element_blank(),
          strip.background    = element_rect(fill = "#DDDDDD", colour = NA),
          strip.text          = element_text(face = "bold", size = 12),
          axis.text.y         = element_text(face = "bold", size = 7),
          axis.title.x        = element_text(face = "bold"),
          plot.title          = element_text(face = "bold", hjust = 0.5),
          panel.clip          = "off"
        )
      # guardar en subcarpeta
      subd <- ont_dirs[[ont]]
      ggsave(
        filename = paste0(t, "_", nm, "_GO_", ont, "_dotplot.png"),
        plot     = p,
        path     = file.path(outdir3, subd),
        width    = 10,
        height   = 6,
        dpi      = 300
      )
      
      # — Aquí añadimos las dos hojas al workbook —
      openxlsx::addWorksheet(wb_go, paste0(ont, "_all"))
      openxlsx::writeData(wb_go, sheet = paste0(ont, "_all"), ego_df)
      openxlsx::addWorksheet(wb_go, paste0(ont, "_top10"))
      openxlsx::writeData(wb_go, sheet = paste0(ont, "_top10"), topn)
    }
    
    # Guardamos el Excel con todos los ontologías
    saveWorkbook(
      wb_go,
      file.path(outdir3, paste0(t, "_", nm, "_GO_enrichment.xlsx")),
      overwrite = TRUE
    )
  }
  message("✔ Enriquecimiento GO completado para ", t)
}

###############################################################################
#   V. ENRIQUECIMIENTO KEGG – PIEL y CEREBRO  (con selección manual de rutas)
###############################################################################
outdir3   <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/3_Enriquecimiento/KEGG"

df_orth  <- readxl::read_excel(orth_file)
mapping  <- unique(df_orth[, c("Dorada_Symbol","Zebra_Entrez")])

conds <- levels(DE_samples[["Piel"]]$coldata$condition)
pairs <- combn(conds, 2, simplify = FALSE)
names(pairs) <- sapply(pairs, \(x) paste0(x[1], "_vs_", x[2]))

# TABLA CON IDS MANUALES  
manual_file <- file.path(outdir3, "manual_KEGG_ids.xlsx")

get_manual_ids <- function(tiss){
  if (!file.exists(manual_file)) return(NULL)
  shs <- getSheetNames(manual_file)
  sheet_to_read <- if (tiss %in% shs) tiss else shs[1]
  
  m <- tryCatch(read.xlsx(manual_file, sheet = sheet_to_read),
                error = \(e) NULL)
  if (is.null(m)) return(NULL)
  
  # Nos aseguramos de que haya 6 columnas exactamente
  colnames(m) <- c("Control_vs_Carragenina",
                   "Control_vs_SST6",
                   "Control_vs_Carragenina_SST6",
                   "Carragenina_vs_SST6",
                   "Carragenina_vs_Carragenina_SST6",
                   "SST6_vs_Carragenina_SST6")
  
  m <- mutate(m, across(everything(), as.character))
  
  m |> mutate(slot = row_number()) |>
    pivot_longer(-slot, names_to = "contrast", values_to = "kegg_id") |>
    filter(nzchar(kegg_id)) |>
    mutate(tissue = tiss)
}

# --- 2) loop tejido/contraste -------------------------------------------------
for (t in names(DE_samples)) {
  
  de_file  <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  man_ids  <- get_manual_ids(t)           # 0 filas ⇢ no hay selección manual
  
  for (nm in names(pairs)) {
    
    # reutilizar o calcular enriquecimiento ---------------------------
    kegg_xlsx <- file.path(outdir3, paste0(t, "_", nm, "_KEGG_enrichment.xlsx"))
    
    if (file.exists(kegg_xlsx)) {
      ekegg_df <- read.xlsx(kegg_xlsx, sheet = "KEGG_all")
    } else {
      res_df <- openxlsx::read.xlsx(de_file, sheet = nm)
      names(res_df)[1] <- "Dorada_Symbol"
      res_df <- res_df[!is.na(res_df$padj), ]
      
      res_ent <- merge(res_df, mapping, by = "Dorada_Symbol") |>
        filter(!is.na(Zebra_Entrez))
      if (nrow(res_ent) == 0) next
      
      geneList <- setNames(res_ent$log2FoldChange,
                           as.character(res_ent$Zebra_Entrez)) |>
        sort(decreasing = TRUE)
      
      ekegg <- gseKEGG(
        geneList     = geneList,
        organism     = "dre",
        minGSSize    = 10,
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )
      ekegg_df <- as.data.frame(ekegg)
      
      # guardar Excel completo (por si no existía)
      wb_tmp <- createWorkbook()
      addWorksheet(wb_tmp, "KEGG_all")
      writeData   (wb_tmp, "KEGG_all", ekegg_df)
      saveWorkbook(wb_tmp, kegg_xlsx, overwrite = TRUE)
    }
    if (nrow(ekegg_df) == 0) next
    
    ekegg_df$Count     <- lengths(strsplit(ekegg_df$core_enrichment, "/"))
    ekegg_df$GeneRatio <- ekegg_df$Count / ekegg_df$setSize
    ekegg_df$Direction <- factor(ifelse(ekegg_df$NES >= 0,
                                        "Activated","Repressed"),
                                 levels = c("Activated","Repressed"))
    
    # selección manual o top-10
    man_sel <- man_ids |> filter(contrast == nm) |> pull(kegg_id)
    
    if (length(man_sel)) {
      topn <- ekegg_df |>
        filter(ID %in% man_sel) |>
        arrange(match(ID, man_sel))       
    } else {
      topn <- ekegg_df |> arrange(p.adjust) |> slice_head(n = 10)
    }
    
    topn <- topn |>
      mutate(Description = str_to_sentence(str_replace_all(Description,"-"," ")))
    
    # dummy-rows para que siempre existan ambos paneles
    if (!"Activated" %in% topn$Direction)
      topn <- bind_rows(topn,
                        data.frame(Description=NA,GeneRatio=NA,Count=NA,
                                   Direction=factor("Activated",
                                                    levels=c("Activated","Repressed"))))
    if (!"Repressed" %in% topn$Direction)
      topn <- bind_rows(topn,
                        data.frame(Description=NA,GeneRatio=NA,Count=NA,
                                   Direction=factor("Repressed",
                                                    levels=c("Activated","Repressed"))))
    
    #  bubble-plot
    p <- ggplot(topn, aes(
      x     = GeneRatio,
      y     = factor(Description, levels = rev(Description)),
      size  = Count,
      fill  = Direction)) +
      geom_point(shape = 21, colour = "black", stroke = 0.5) +
      scale_fill_manual(
        name   = "Direction",
        values = c(Activated = "#fb8072", Repressed = "#80b1d3"),
        breaks = c("Activated","Repressed"), drop = FALSE) +
      guides(fill = guide_legend(
        override.aes = list(size = 5, shape = 21, colour = "black",
                            fill = c("#fb8072","#80b1d3")),
        order = 1, title.position = "top")) +
      scale_size_continuous(range = c(3,8)) +
      facet_grid(. ~ Direction, scales = "free_x", drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.06,0.06))) +
      coord_cartesian(clip = "off") +
      scale_y_discrete(na.translate = FALSE) +
      labs(title = paste(t, nm, "KEGG", sep = " | "),
           x = "Gene ratio", y = NULL) +
      theme_bw(base_size = 14) +
      theme(panel.spacing = unit(1,"lines"),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "#DDDDDD", colour = NA),
            strip.text       = element_text(face = "bold", size = 12),
            axis.text.y      = element_text(face = "bold", size = 7),
            axis.title.x     = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5),
            panel.clip       = "off")
    
    ggsave(
      filename = paste0(t, "_", nm, "_GSEA_KEGG_dotplot.png"),
      plot     = p,
      path     = outdir3,
      width    = 10, height = 6, dpi = 300)
    
    # añade hoja Top-10 si aún no existe -----------------------------
    if (!file.exists(kegg_xlsx) ||
        !"KEGG_top10" %in% getSheetNames(kegg_xlsx)) {
      wb_kegg <- loadWorkbook(kegg_xlsx)
      addWorksheet(wb_kegg, "KEGG_top10")
      writeData   (wb_kegg, "KEGG_top10", topn)
      saveWorkbook(wb_kegg, kegg_xlsx, overwrite = TRUE)
    }
  }
  message("✔ GSEA KEGG completado para ", t)
}

###############################################################################
# VI. ENRIQUECIMIENTO Reactome – PIEL y CEREBRO (con selección manual de rutas)
###############################################################################
outdir4 <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/3_Enriquecimiento/Reactome"

# tabla de IDS manuales
manual_file <- file.path(outdir4, "manual_reactome_ids.xlsx")

get_manual_ids <- function(tiss){
  if (!file.exists(manual_file)) return(NULL)
  
  shs <- getSheetNames(manual_file)
  sheet_to_read <- if (tiss %in% shs) tiss else shs[1]
  
  m <- tryCatch(read.xlsx(manual_file, sheet = sheet_to_read),
                error = \(e) NULL)
  if (is.null(m)) return(NULL)
  
  ## nos aseguramos de que haya exactamente 6 columnas
  colnames(m) <- c("Control_vs_Carragenina",
                   "Control_vs_SST6",
                   "Control_vs_Carragenina_SST6",
                   "Carragenina_vs_SST6",
                   "Carragenina_vs_Carragenina_SST6",
                   "SST6_vs_Carragenina_SST6")
  
  ## convertimos TODO a character para evitar mezclas de tipos
  m <- dplyr::mutate(m, dplyr::across(everything(), as.character))
  
  m |>
    mutate(slot = row_number()) |>
    tidyr::pivot_longer(-slot,
                        names_to  = "contrast",
                        values_to = "reactome_id") |>
    filter(nzchar(reactome_id)) |>
    mutate(tissue = tiss)
}

# bucle principal 
for (t in names(DE_samples)) {
  
  de_file <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  man_ids <- get_manual_ids(t)      # 0 filas ⇢ sin selección manual
  
  for (nm in names(pairs)) {
    
    ## ENRIQUECIMIENTO
    react_xlsx <- file.path(outdir4, paste0(t, "_", nm, "_Reactome_enrichment.xlsx"))
    
    if (file.exists(react_xlsx)) {
      # ▸ ya existe → lee hoja Reactome_all
      ereact_df <- read.xlsx(react_xlsx, sheet = "Reactome_all")
    } else {
      # ▸ no existe → se calcula como siempre (tu código original)
      res_df <- openxlsx::read.xlsx(de_file, sheet = nm)
      names(res_df)[1] <- "Dorada_Symbol"
      res_df <- res_df[!is.na(res_df$padj), ]
      res_ent <- merge(res_df, mapping, by = "Dorada_Symbol") |>
        filter(!is.na(Zebra_Entrez))
      if (nrow(res_ent) == 0) next
      
      geneList <- setNames(res_ent$log2FoldChange,
                           as.character(res_ent$Zebra_Entrez)) |>
        sort(decreasing = TRUE)
      
      ereact <- ReactomePA::gsePathway(
        geneList     = geneList,
        organism     = "zebrafish",
        minGSSize    = 10,
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )
      ereact_df <- as.data.frame(ereact) |>
        mutate(Direction = ifelse(NES >= 0,"Activated","Repressed"))
      
      # Excel completo (por si no existía)
      wb_tmp <- openxlsx::createWorkbook()
      addWorksheet(wb_tmp, "Reactome_all")
      writeData   (wb_tmp, "Reactome_all", ereact_df)
      saveWorkbook(wb_tmp, react_xlsx, overwrite = TRUE)
    }
    if (nrow(ereact_df) == 0) next
    
    ereact_df$Count     <- lengths(strsplit(ereact_df$core_enrichment, "/"))
    ereact_df$GeneRatio <- ereact_df$Count / ereact_df$setSize
    
    ## Selección manual del TOP-10
    man_sel <- man_ids |>
      filter(contrast == nm) |>
      pull(reactome_id)
    
    if (length(man_sel)) {
      topn <- ereact_df |>
        filter(ID %in% man_sel) |>
        arrange(match(ID, man_sel))
    } else {
      topn <- ereact_df |>
        arrange(p.adjust) |>
        slice_head(n = 10)
    }
    
    topn <- topn |>
      mutate(
        Description = str_to_sentence(str_replace_all(Description,"-"," ")),
        Direction   = factor(Direction, levels = c("Activated","Repressed"))
      )
    
    # dummy-rows
    if (!"Activated" %in% topn$Direction)
      topn <- bind_rows(topn, data.frame(Description = NA, GeneRatio = NA,
                                         Count = NA,
                                         Direction = factor("Activated",
                                                            levels = c("Activated","Repressed"))))
    if (!"Repressed" %in% topn$Direction)
      topn <- bind_rows(topn, data.frame(Description = NA, GeneRatio = NA,
                                         Count = NA,
                                         Direction = factor("Repressed",
                                                            levels = c("Activated","Repressed"))))
    
    ##  Bubble-plot
    p <- ggplot(topn, aes(
      x     = GeneRatio,
      y     = factor(Description, levels = rev(Description)),
      size  = Count,
      fill  = Direction
    )) +
      geom_point(shape = 21, colour = "black", stroke = 0.5) +
      scale_fill_manual(
        name   = "Direction",
        values = c(Activated = "#fb8072", Repressed = "#80b1d3"),
        breaks = c("Activated","Repressed"), drop = FALSE
      ) +
      guides(fill = guide_legend(
        override.aes = list(size = 5, shape = 21, colour = "black",
                            fill = c("#fb8072","#80b1d3")),
        order = 1, title.position = "top"
      )) +
      scale_size_continuous(range = c(3,8)) +
      facet_grid(. ~ Direction, scales = "free_x", drop = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(0.06,0.06))) +
      coord_cartesian(clip = "off") +
      scale_y_discrete(na.translate = FALSE) +
      labs(
        title = paste(t, nm, "Reactome", sep = " | "),
        x     = "Gene ratio", y = NULL
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.spacing    = unit(1,"lines"),
        panel.background = element_rect(fill = "white"),
        panel.border     = element_rect(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#DDDDDD", colour = NA),
        strip.text       = element_text(face = "bold", size = 12),
        axis.text.y      = element_text(face = "bold", size = 7),
        axis.title.x     = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        panel.clip       = "off"
      )
    
    ggsave(
      filename = paste0(t, "_", nm, "_GSEA_Reactome_dotplot.png"),
      plot     = p,
      path     = outdir4,
      width    = 10, height = 6, dpi = 300
    )
    
    ## añade hoja top-10 (solo si no existía) 
    if (!file.exists(react_xlsx) || !"Reactome_top10" %in%
        getSheetNames(react_xlsx)) {
      wb_react <- loadWorkbook(react_xlsx)
      addWorksheet(wb_react, "Reactome_top10")
      writeData   (wb_react, "Reactome_top10", topn)
      saveWorkbook(wb_react, react_xlsx, overwrite = TRUE)
    }
  }
  message("✔ GSEA Reactome completado para ", t)
}

##############################################################################
# Paneles resumen KEGG + Reactome
###############################################################################
# Paneles con 10 términos cada uno
manual_panels <- tribble(
  ~tissue,   ~contrast,                          ~terms,
  "Cerebro", "Control_vs_Carragenina", 
  c(
    "Tak1 dependent ikk and nf kappa b activation",
    "Neuron development",
    "Cell cell signaling",
    "Sumoylation of ubiquitinylation proteins",
    "Regulation of pten stability and activity",
    "Pten regulation",
    "Primary bile acid biosynthesis",
    "Neurotransmitter release cycle",
    "Chromosome condensation",
    "Apoptotic execution phase"
  ),
  "Piel",    "Control_vs_Carragenina",
  c(
    "Antigen processing: ubiquitination & proteasome degradation",
    "Positive regulation of t cell activation",
    "Chromosome condensation",
    "Metabolic pathways",
    "Neurogenesis",
    "Keratinization",
    "Cytokine cytokine receptor interaction",
    "Apoptosis",
    "Pten regulation",
    "Innate immune system"
  ),
  "Cerebro", "Control_vs_SST6",
  c(
    "Neurotransmitter transport",
    "Pten regulation",
    "Regulation of wnt signaling pathway",
    "Regulation of rna biosynthetic process",
    "Sensory organ development",
    "Neurotransmitter release cycle",
    "Cell cell signaling",
    "Antigen processing: ubiquitination & proteasome degradation",
    "Brain development",
    "Keratinization"
  ),
  "Piel",    "Control_vs_SST6",
  c(
    "Antigen processing and presentation",
    "Keratinization",
    "Circadian rhythm",
    "Response to stress",
    "Regulation of biosynthetic process",
    "Lipid metabolic process",
    "Apoptotic process",
    "Innate immune system",
    "Neutrophil degranulation",
    "Metabolic pathways"
  ),
  # NUEVAS comparaciones Control_vs_Carragenina_SST6:
  "Cerebro", "Control_vs_Carragenina_SST6",
  c(
    "Programmed cell death",
    "Pten Regulation",
    "Tak1 dependent ikk and nf kappa b activation",
    "Interleukin 1 signaling",
    "Cell cell communication",
    "Sumoylation of ubiquitinylation proteins",
    "Neurotransmitter release cycle",
    "Response to growth factor",
    "Primary bile acid biosynthesis",
    "Regulation of nf kappa b signaling"
  ),
  "Piel",    "Control_vs_Carragenina_SST6",
  c(
    "Adaptive immune system",
    "Neuronal system",
    "Apoptosis",
    "Pten regulation",
    "Response to stress",
    "Lymphocyte activation",
    "Epidermis development",
    "Signaling by wnt",
    "Ppar signaling pathway",
    "Cytoskeleton in muscle cells"
  )
)

# Helper para cargar y normalizar cada GSEA
load_gsea <- function(tissue, contrast, db, enr_root) {
  base_dir <- switch(db,
                     "GO-BP"    = enr_root,
                     "KEGG"     = file.path(enr_root, "KEGG"),
                     "Reactome" = file.path(enr_root, "Reactome"))
  suffix   <- switch(db,
                     "GO-BP"    = "GO",
                     "KEGG"     = "KEGG",
                     "Reactome" = "Reactome")
  file <- file.path(base_dir,
                    paste0(tissue, "_", contrast, "_", suffix, "_enrichment.xlsx"))
  sheet <- switch(db,
                  "GO-BP"    = "BP_all",
                  "KEGG"     = "KEGG_all",
                  "Reactome" = "Reactome_all")
  if (!file.exists(file)) return(tibble())
  
  df0 <- readxl::read_xlsx(file, sheet = sheet)
  
  # Count: prefer core_enrichment, luego columna Count, sino 0
  if ("core_enrichment" %in% names(df0)) {
    df0$Count <- lengths(strsplit(df0$core_enrichment, "/"))
  } else if ("Count" %in% names(df0)) {
    df0$Count <- as.numeric(df0$Count)
  } else {
    df0$Count <- 0
  }
  
  # GeneRatio: prefer columna GeneRatio, luego setSize
  if ("GeneRatio" %in% names(df0)) {
    df0$GeneRatio <- as.numeric(df0$GeneRatio)
  } else if ("setSize" %in% names(df0)) {
    df0$GeneRatio <- df0$Count / df0$setSize
  } else {
    df0$GeneRatio <- NA_real_
  }
  
  # Dirección según NES
  df0$Direction <- factor(
    ifelse(df0$NES >= 0, "Activated", "Repressed"),
    levels = c("Activated","Repressed")
  )
  df0$Source <- db
  df0$Description <- stringr::str_to_sentence(
    str_replace_all(df0$Description, "-", " ")
  )
  
  df0 %>%
    dplyr::select(Source, Description, GeneRatio, Count, Direction)
}

# Construcción de los 6 gráficos
enr_root <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/3_Enriquecimiento"
shapes   <- c("GO-BP"=21, "KEGG"=22, "Reactome"=24)

plots <- manual_panels %>%
  pmap(function(tissue, contrast, terms) {
    df_all <- map_dfr(c("GO-BP","KEGG","Reactome"),
                      ~ load_gsea(tissue, contrast, .x, enr_root))
    df <- df_all %>% filter(Description %in% terms)
    
    grid <- expand_grid(
      Direction   = factor(c("Activated","Repressed"),
                           levels = c("Activated","Repressed")),
      Description = factor(terms, levels = rev(terms)),
      Source      = factor(c("GO-BP","KEGG","Reactome"),
                           levels = c("GO-BP","KEGG","Reactome"))
    )
    
    p <- ggplot(
      grid %>% left_join(df, by = c("Source","Direction","Description")),
      aes(x = GeneRatio, y = Description,
          size  = Count,
          shape = Source,
          fill  = Direction)
    ) +
      geom_point(
        colour = "black",         
        stroke = 0.4,         
        na.rm  = TRUE
      ) +
      
      # shape según base de datos
      scale_shape_manual(
        values = shapes,
        name   = "Database"
      ) +
      
      # fill según activated/repressed
      scale_fill_manual(
        values = c(Activated = "#fb8072", Repressed = "#80b1d3"),
        name   = "Direction"
      ) +
      
      # size → círculos blancos en la leyenda
      scale_size_continuous(
        range = c(2, 6),
        name  = "Count",
        guide = guide_legend(
          override.aes = list(
            shape = 21,
            fill  = "white",
            colour= "black"
          ),
          title.position = "top",
          order = 3
        )
      ) +
      
      facet_grid(
        . ~ Direction,
        scales = "fixed",
        space  = "fixed"
      ) +
      
      labs(
        title = paste0(tissue, " – ", gsub("_", " ", contrast)),
        x     = "Gene ratio",
        y     = NULL
      ) +
      
      scale_x_continuous(
        expand = expansion(mult = c(0.2, 0.2))
      ) +
      scale_y_discrete(
        expand = expansion(mult = c(0.1, 0.1))
      ) +
      
      # IMPORTANTÍSIMO para que no se recorten puntos ni cuadros
      coord_cartesian(clip = "off") +
      
      theme_bw(base_size = 12) +
      theme(
        # márgenes extra por la derecha
        plot.margin = margin(t = 5, r = 100, b = 5, l = 5, unit = "pt"),
        
        # ubicación fija de la leyenda
        legend.position      = c(1.02, 0.50),
        legend.justification = c("left", "center"),
        legend.direction     = "vertical",
        
        legend.box.margin    = margin(0, 0, 0, 0),
        legend.margin        = margin(5, 5, 5, 5),
        legend.key.size      = unit(0.8, "lines"),
        legend.title         = element_text(size = 10, face = "bold"),
        legend.text          = element_text(size = 8),
        
        # paneles y títulos
        plot.title       = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "#DDDDDD", colour = NA),
        strip.text       = element_text(face = "bold", size = 10),
        axis.text.y      = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      
      # Para garantizar orden y consistencia en la leyenda de shapes/fill
      guides(
        shape = guide_legend(
          title.position = "top",
          order = 1
        ),
        fill  = guide_legend(
          title.position = "top",
          override.aes  = list(shape = 21, colour = "black"),
          order = 2
        )
      )
    
    list(key = paste(tissue, contrast, sep = "___"), plot = p)
  })

# Finalmente se guardan
outdir <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/4_GraficosYAnalisisFinales/Main_figures"
walk(plots, ~ ggsave(
  filename = file.path(outdir, paste0(.x$key, ".png")),
  plot     = .x$plot,
  width    = 6,
  height   = 4,
  dpi      = 300
))

###############################################################################
#+ Pathview de KEGG (ruta fijada por contraste)  +  Excel con genes coloreados
###############################################################################
pv_dir <- file.path(outdir5, "Pathview")
dir.create(pv_dir, showWarnings = FALSE)

forced <- data.frame(
  tissue   = rep(c("Piel", "Cerebro"), each = 6),
  contrast = rep(c("Control_vs_Carragenina",
                   "Control_vs_SST6",
                   "Control_vs_Carragenina_SST6",
                   "Carragenina_vs_SST6",
                   "Carragenina_vs_Carragenina_SST6",
                   "SST6_vs_Carragenina_SST6"), times = 2),
  kegg_id  = c(              # ← pon aquí las vías fijas que quieras
    "04060", "01100", "04820", "04216", "04820", "04216",  # Piel
    "00120", "", "00120", "00860", "00860", ""  # Cerebro
  ),
  stringsAsFactors = FALSE
)

for (t in names(DE_samples)) {
  
  xls_DE <- file.path(outdir2, paste0(t, "_DESeq2_all_comparisons.xlsx"))
  shs    <- getSheetNames(xls_DE)[grepl("_vs_", getSheetNames(xls_DE))]
  
  for (sh in shs) {
    
    ## rutas forzadas
    row_f <- forced[forced$tissue == t & forced$contrast == sh, , drop = FALSE]
    pw_code <- if (nrow(row_f) == 1 && nzchar(row_f$kegg_id)) row_f$kegg_id else NA
    
    ## Si no hay ruta forzada, toma la primera ≠01100 de KEGG_top10
    if (is.na(pw_code)) {
      kegg_xlsx <- file.path(kegg_dir, paste0(t, "_", sh, "_KEGG_enrichment.xlsx"))
      if (!file.exists(kegg_xlsx)) next                # no hay enriquecimiento
      top10 <- tryCatch(read.xlsx(kegg_xlsx, sheet = "KEGG_top10"),
                        error = \(e) NULL)
      if (is.null(top10) || nrow(top10) == 0) next
      valid <- top10$ID[!is.na(top10$ID) & top10$ID != "" & top10$ID != "dre01100"]
      if (!length(valid)) next
      pw_code <- sub("^dre", "", valid[1])
    }
    
    ## Vector ENTREZ → log2FC
    df <- read.xlsx(xls_DE, sheet = sh)
    gL <- setNames(df$log2FoldChange, df$gene_id)
    names(gL) <- mapping$Zebra_Entrez[match(names(gL), mapping$Dorada_Symbol)]
    gL <- gL[!is.na(names(gL))]
    if (!length(gL)) next
    
    message(sprintf("→ Pathview %-8s | %-35s | dre%s", t, sh, pw_code))
    
    out_xlsx <- file.path(outdir5,
                          sprintf("%s_%s_Pathview_gene_values.xlsx", t, sh))
    
    ## Pathview
    oldwd <- setwd(pv_dir)
    res <- tryCatch(
      pathview(gene.data  = gL,
               cpd.data   = NULL,
               pathway.id = pw_code,
               species    = "dre",
               out.suffix = paste(t, sh, sep = "_"),
               limit      = list(gene = 3, cpd = 1)),
      error = \(e) {message("  ⚠  Pathview falló: ", e$message); NULL}
    )
    setwd(oldwd)
    
    ## Excel con genes coloreados
    if (is.null(res) || is.null(res$plot.data.gene)) {
      write.xlsx(data.frame(), out_xlsx, overwrite = TRUE)
    } else {
      write.xlsx(
        data.frame(
          ENTREZ = rownames(res$plot.data.gene),
          log2FC = res$plot.data.gene$mol.data
        ),
        out_xlsx, overwrite = TRUE, rowNames = FALSE
      )
    }
  }
}

























###############################################################################
###### Para hacer la tabla de anotación para analizar los sesgos
###############################################################################
library(GenomicFeatures)
library(Rsamtools)
library(Biostrings)

gtf    <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/Annotations.gtf"
fasta  <- "/home/alumno02/TFM/AnalisisDeRNAseq/data/Genome_and_annotations/GCF_900880675.1_fSpaAur1.1_genomic.fna"
out.tsv <- "/home/alumno02/TFM/AnalisisDeRNAseq/annotation_dorada.tsv"

# 1. Construye un TxDb a partir del GTF
txdb <- makeTxDbFromGFF(gtf, format = "gtf")

# 2. Rango exónico por gen y longitud (sin solapamientos)
exByGene <- exonsBy(txdb, by = "gene")   # lista GRanges por gen
exRed    <- reduce(exByGene)             # exones fusionados dentro de cada gen

## Longitud génica = suma de anchos de los exones fusionados
gene_len <- sum(width(exRed))            # named integer vector


# 3. GC por gen: recorrer la lista (cada elemento es GRanges)

fa <- FaFile(fasta); open(fa)            # asegúrate de tener el .fai

gc_content <- sapply(names(exRed), function(g) {
  seqs <- getSeq(fa, exRed[[g]])         # exRed[[g]] es GRanges → OK
  # Cuenta G y C en todas las exones y divide por la longitud total
  sum(letterFrequency(seqs, c("G","C"))) / sum(width(seqs))
})

# 4. Empaqueta y guarda
annotation <- data.frame(
  gene_id = names(gene_len),
  length  = as.numeric(gene_len),
  gc      = round(gc_content, 4),
  row.names = NULL
)

write.table(annotation,
            file = out.tsv,
            quote = FALSE, sep = "\t", row.names = FALSE)

cat("Tabla de anotación creada en:\n", out.tsv, "\n")

## Descarta genes con NA/NaN en length o gc
annotation <- na.omit(annotation)

cat(nrow(annotation), "genes conservados tras eliminar NA/NaN\n")

write.table(annotation,
            file = out.tsv,
            quote = FALSE, sep = "\t", row.names = FALSE)

cat("Tabla de anotación creada en:\n", out.tsv, "\n")








###############################################################################
###### Para generar el archivo de los Ortologos Dorada y Pez cebra
###############################################################################
# 1) Carga de librerías
library(readxl)
library(biomaRt)
library(dplyr)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(openxlsx)

# 2) Leer la tabla desde Excel
Ortologos_pez_zebra <- read_excel(
  "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/1_DEG_Dorada.xlsx"
)

# 3) Conéctarse a Ensembl (zebrafish)
mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# 4) Preparar los ENSG (quita NAs y duplicados)
ids <- unique(na.omit(Ortologos_pez_zebra$ortholog_ensg))

# 5) Consultar terminos GO
go_annot <- getBM(
  attributes = c(
    "ensembl_gene_id",  # ENSDARG...
    "go_id",            # GO accession
    "name_1006",        # GO term name
    "namespace_1003"    # GO aspect (BP/MF/CC)
  ),
  filters = "ensembl_gene_id",
  values  = ids,
  mart    = mart
)

# 6) Fusionar de nuevo con la tabla original
Ortologos_pez_zebra_con_GO <- merge(
  Ortologos_pez_zebra,
  go_annot,
  by.x  = "ortholog_ensg",
  by.y  = "ensembl_gene_id",
  all.x = TRUE
)

# Renombrar las columnas GO para mayor claridad
Ortologos_pez_zebra_con_GO <- Ortologos_pez_zebra_con_GO %>%
  rename_with(~ recode(.x,
                       initial_alias    = "Dorada_Symbol",
                       initial_ensg     = "Dorada_Ensembl",
                       `o#`             = "Ortholog_Order",
                       ortholog_name    = "Zebra_Symbol",
                       ortholog_ensg    = "Zebra_Ensembl",
                       description      = "Description",
                       go_id            = "GO_ID",
                       name_1006        = "GO_Term",
                       namespace_1003   = "GO_Category"
  ))


# 7) Renombrar y ordenar
Ortologos_pez_zebra_con_GO <- 
  Ortologos_pez_zebra_con_GO[
    order(Ortologos_pez_zebra_con_GO$Ortholog_Order),
  ]

# 8) Echar un vistazo
head(Ortologos_pez_zebra_con_GO)

# 9) Guardar en excel
openxlsx::write.xlsx(
  Ortologos_pez_zebra_con_GO,
  "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/2_Ortologos_pez_zebra_con_GO.xlsx",
  rowNames = FALSE
)

# 10) Añadir columna ENTREZID de Zebrafish para facilitar el downstream
# 10.1) Obtenemos el mapeo Ensembl→Entrez y quitamos duplicados
mapping <- AnnotationDbi::select(
  org.Dr.eg.db,
  keys    = unique(Ortologos_pez_zebra_con_GO$Zebra_Ensembl),
  keytype = "ENSEMBL",
  columns = "ENTREZID"
)

# Nos quedamos con el primer ENTREZID si hay duplicados
mapping2 <- mapping[!duplicated(mapping$ENSEMBL), ]

# 10.2) Unimos a tu tabla original, usando base R para evitar sorpresas
Ortologos_pez_zebra_con_GO <-
  merge(
    Ortologos_pez_zebra_con_GO,
    mapping2,
    by.x  = "Zebra_Ensembl",
    by.y  = "ENSEMBL",
    all.x = TRUE
  )

# 10.3) Renombramos la columna que ahora sí existe
names(Ortologos_pez_zebra_con_GO)[
  names(Ortologos_pez_zebra_con_GO) == "ENTREZID"
] <- "Zebra_Entrez"

# Y comprobamos
head(Ortologos_pez_zebra_con_GO)


# 10.4) Guardar el resultado final en Excel
openxlsx::write.xlsx(
  Ortologos_pez_zebra_con_GO,
  file      = "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/3_Ortologos_Zebra_con_Entrez_GO.xlsx",
  rowNames  = FALSE,
  overwrite = TRUE
)

# 11) Creamos el objeto final con los nombres actualizados
DEG_DoradaZebra_GO <- Ortologos_pez_zebra_con_GO[
  order(Ortologos_pez_zebra_con_GO$Ortholog_Order),
  c(
    "Dorada_Symbol",    # símbolo en dorada
    "Dorada_Ensembl",   # ENSG dorada
    "Ortholog_Order",   # tu orden original
    "Zebra_Symbol",     # símbolo ortólogo pez cebra
    "Zebra_Ensembl",    # ENSDARG ortólogo
    "Zebra_Entrez",     # ENTREZID pez cebra
    "Description",      # descripción
    "GO_ID",            # GO term
    "GO_Term",          # nombre GO
    "GO_Category"       # categoría GO (BP/MF/CC)
  )
]

# 12) Guardamos el resultado final
openxlsx::write.xlsx(
  DEG_DoradaZebra_GO,
  file      = "/home/alumno02/TFM/AnalisisDeRNAseq/data/6_DEG/0_Generacion_Archivo_Ortologos_Dorada/4_DEG_DoradaZebra_GO.xlsx",
  rowNames  = FALSE,
  overwrite = TRUE
)
