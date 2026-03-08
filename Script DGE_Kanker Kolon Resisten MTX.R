#CAPSTONE PROJECT OMICS LITE | NATALIA RARA 423
#Platform: Affymetrix Human Genome U133 Plus 2.0 Array
#Dataset: GSE11440 (Human HT29 Colon Cancer, Sensitive vs Resistance to Methotrexate (MTX))
#Tujuan: Mengetahui gen yang diekspresikan berbeda di sel yang resisten MTX

#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen
#antara sel kanker yang resisten MTX dan sel kanker sensitive MTX
#Analisis DEG dilakukan dengan pendekatan LIMMA

#Instal platorm hgu133plus2
BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)

#MEMANGGIL LIBRARY UNTUK MENGAKTIFKAN PACKAGE
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)

#PENGAMBILAN DATA DARI GEO
gset <- getGEO("GSE11440", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PREPOSESING DATA EKSPRESI
ex <- exprs(gset)
dim(ex)

#normalisasi data
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
dim(ex) 
qx

LogTransform <- (qx[5] > 100) | (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
dim(ex) 
qx

#MEMBERI DEFINISI KELOMPOK SAMPEL
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)

#MEMBUAT DESAIN MATRIX
design <- model.matrix(~0 + gset$group)
#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)

#Menentukan perbandingan biologis
grup_resisten <- nama_grup[1]
grup_sensitive <- nama_grup[2]
contrast_formula <- paste(grup_resisten, "-", grup_sensitive)
print(paste("Kontras yang dianalisis:", contrast_formula))

#ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

#Mengambil hasil deg
deg <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(deg)

#memfilter deg signifikan
deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
head(deg_sig)
nrow(deg_sig)

#ANOTASI GEN
probe_ids <- rownames(deg)

gene_annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#menerapkan hasil anotasi ke deg
deg$PROBEID <- rownames(deg)
deg <- merge(
  deg,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(deg[, c("PROBEID", "SYMBOL", "GENENAME")])

#MEMBUAT BOXPLOT DISTRIBUSI NILAI EKSPRESI
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.7
)

#DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

#UMAP (VISUALISASI DIMENSI RENDAH)
#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#mengecek apakah di umap_input ada NA
any(is.na(umap_input))

#kalau TRUE, lanjut imputasi pakai 0
umap_input[is.na(umap_input)] <- 0

#Jalankan UMAP
umap_result <- umap(umap_input, n_neighbors = 2)  # < 6


#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#VISUALISASI VOLCANO PLOT
volcano_data <- data.frame(
  logFC = deg$logFC,
  adj.P.Val = deg$adj.P.Val,
  Gene = deg$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Kolon Resisten vs Sensitive")

#VISUALISASI HEATMAP
mat_heatmap <- ex[deg$PROBEID, ]
#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(deg$SYMBOL) | deg$SYMBOL == "",
  deg$PROBEID, # jika SYMBOL kosong → probe ID
  deg$SYMBOL # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap
pheatmap(
  mat_heatmap,
  scale = "row", # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE, # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Significant Differentially Expressed Genes"
)

#MENYIMPAN HASIL
# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(deg, "Hasil_GSE11440_DEG.csv")
message("Analisis selesai. File hasil telah disimpan.")

#ANALISIS ENRICHMENT GO (FUNGSI BIOLOGIS) & KEGG
#TAHAPAN ENRICHMENT ANALYSIS
genes <- unique(deg$SYMBOL)
head(genes)

#Instalasi dan loading library
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

#koknversi ke ENTREZ
#KEGG butuh bentuk ENTREZ ID (numerik) bukan PROBEID
gene_entrez <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#cek hasil konversi
head(gene_entrez)

#Analisis GO Biological Process
ego_bp <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

#melihat hasil GO
head(ego_bp)

#Analisis KEGG
ekegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa", # human
  pvalueCutoff = 0.05
)

#melihat hasil KEGG
head(ekegg)

#Visualisasi GO dan KEGG dalam bentuk dotplot
#visualisasi GO
dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process")

barplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment")

#visualisasi KEGG
dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway")

#menyimpan file ke Excel
write.csv(as.data.frame(ego_bp), "GO_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_results.csv")