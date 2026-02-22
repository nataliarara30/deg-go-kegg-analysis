**ANALISIS DIFFERENTIALLY EXPRESSED GENES (DEG) KANKER KOLOREKTAL (PAIRED DESIGN), GO & KEGG**

1. **PENDAHULUAN**

   Kanker kolorektal merupakan salah satu jenis kanker dengan tingkat mortalitas yang tinggi di seluruh dunia. Perkembangan teknologi genomik memungkinkan analisis ekspresi gen secara simultan dalam skala besar menggunakan data microarray. Analisis DEG dapat digunakan unutk mengidentifikasi gen-gen yang berperan dalam perkembangan dan progresi kanker.

   Project ini bertujuan mengidentifikasi gen-gen yang mengalami perubahan ekspresi pada kanker kolorektal adenocarcinoma menggunakan database GSE110223 dari GEO, serta menganalisis fungsi biologis dan jalur molekuler yang terlibat melalui analisis Gene Ontology (GO) dan Kyoto Encyclopedia of Genes and Genomes (KEGG)

1. **METODE**
1. Pengambilan Data

   Data ekspresi gen menggunakan dataset GSE110223 dari database GEO menggunakan package GEOquery pada software R 4.5.2. Dataset yang digunakan merupakan data microarray berbasis platform Affymetrix yang membandingkan jaringan kanker kolorektal dan jaringan normal.

1. Preparasi Data

   Data ekspresi diekstraksi menggunakan fungsi exprs(). Selanjutnya dilakukan pemeriksaan distribusi data menggunakan boxplot untuk memastikan kualitas data. Transformasi log2 dilakukan apabila diperlukan berdasarkan distribusi nilai ekspresi.

1. Analisis DEG

   Analisis DEG dilakukan menggunakan paket limma. Matriks desain (design matrix) dibangun berdasarkan kelompok sampel tumor (tumor dan normal). Model linear kemudian difitting menggunakan fungsi lmfit() dan eBayes().

   Gen dianggap signifikan apabila memenuhi kriteria:

- adjusted p-value < 0,05
- |logFC| > 1


1. Anotasi Gen

   ID probe Affymetrix dianotasi menggunakan package anotasi yang sesuai. ID probe dikonversi menjadi symbol dan nama gen.

1. Visualisasi Data

   Visualisasi hasil dilakukan menggunakan:

- Volcano plot untuk menampilkan distribusi DEG
- Heatmap untuk menampilkan pola ekspresi 50 gen paling signifikan
- Dotplot GO dan KEGG untuk analisis fungsional

1. Analisis GO dan KEGG

   Analisis Gene Ontology dan KEGG dilakukan menggunakan paket clusterProfiler. Simbol gen dikonversi menjadi Entrez ID menggunakan org.Hs.eg.db. Enrichment analysis dilakukan dengan batas signifikansi p-adjust < 0,05.

1. **HASIL DAN INTERPRETASI**
1. Volcano Plot

   Volcano plot (Gambar 1) menunjukkan distribusi gen yang mengalami perubahan ekspresi antara janringa kanker dan normal. gen dengan nilai logFC > 1 dan adjusted p-value < 0,05 dikategorikan sebagai upregulated genes (merah), sedangkan  gen logFC < -1 dan adjusted p-value < 0,05 dikategorikan sebagai downregulated genes (biru).

   Hasil analisis menunjukkan bahwa terdapat 733 gen yang mengalami perubahan ekspresi signifikan baik upregulated maupun downregulated. 

   ![volcano plot](Aspose.Words.9b11e175-1a97-4fb6-9c1c-3eb785629a2d.001.png)

   Gambar 1. Volcano Plot.

1. 50 Gen Signifikan

   Heatmap (Gambar 2) menampilkan pola ekspresi 50 gen paling signifikan berdasarkan nilai adjusted p-value. Warna merah menunjukkan tingkat ekspresi gen yang tinggi, sedangkan warna biru menunjukkan tingkat ekspresi gen yang rendah. Hasil clustering menunjukkan bahwa sampel kanker dan sampel normal membentuk dua kelompok utama yang terpisah secara jelas. Sampel kanker kolorektal menunjukkan pola ekspresi yang berbeda secara konsisten dibandingkan dengan jaringan normal. Beberapa gen terlihat mengalami peningkatan ekspresi pada kelompok kanker, sementara gen lainnya mengalami penurunan ekspresi. Pemisahan ini menunjukkan bahwa gen-gen tertentu dapat digunakan sebagai biomarker molekuler adanya kanker kolorektal.

   ![heatmap](Aspose.Words.9b11e175-1a97-4fb6-9c1c-3eb785629a2d.002.png)

   Gambar 2. Heatmap.

1. Analisis GO Fungsi Biologis

   Hasil analisis GO fungsi biologis (Gambar 3) menunjukkan bahwa gen-gen yang terekspresi berbeda banyak terlibat dalam proses biologis berikut seperti respons terhadap stimulus xenobiotik, proses metabolisme steroid, respons terhadap level oksigen, regulasi proses metabolisme molekul kecil, respons terhadap substansi toksik, dan detoksifikasi. Dominasi proses respon terhadap xenobiotik dan toksin menunjukkan bahwa sel  kanker mengalami perubahan sistem detoksifikasi dan metabolisme zat asing. Hal ini berkaitan dengan adaptasi sel kanker terhadap lingkungan yang stres dan kondisi metabolik yang tidak stabil. Selain itu, keterlibatan proses respons terhadap hipoksia dan kadar oksigen rendah menunjukkan  adanya adaptasi terhadap kondisi microenvironment tumor yang sering mengalami kekurangan oksigen. Proses metabolisme hormon dan steroid juga mengindikasikan adanya perubahan regulasi hormonal yang berperan dalam progresi kanker. 

   ![GO dotplot](Aspose.Words.9b11e175-1a97-4fb6-9c1c-3eb785629a2d.003.png)

   Gambar 3. GO fungsi biologis dotplot.

1. Analisis KEGG Pathway

   Hasil analisis KEGG (Gambar 4) menunjukkan  bahwa gen-gen yang terekspresi berbeda terlibat dalam beberapa pathway metabolik dan regulasi sel seperti karsinogenesis kimia dan aktivasi reseptor, siklus sel, sekresi empedu, biosintesis kofaktor, biosintesis hormon steroid, metabolisme retinol, metabolisme obat (sitokrom P450), dan metabolisme xenobiotik. 

   ![KEGG pathway dotplot](Aspose.Words.9b11e175-1a97-4fb6-9c1c-3eb785629a2d.004.png)

   Gambar 4. KEGG Pathway dotplot.

   Pathway siklus sel menunjukkan  bahwa gen-gen signifikan berperan dalam regulasi pembelahan sel yang merupakan karakteristik utama sel kanker dan berkaitan dengan proliferasi yang tidak terkendali. Jalur metabolisme xenobiotik dan sitokrom P450 menunjukkan perubahan kemampuan sel kanker dalam memetabolisme obat dan zat asing yang berpotensi berkontribusi terhadap resistensi terapi. Selain itu, keterlibatan jalur karsinogenesis kimia mengindikasikan bahwa gen-gen tersebut berperan dalam proses aktivasi karsinogen yagn dapat memicu mutasi dan perkembangan kanker.

1. **KESIMPULAN**
- Analisis DEG menggunakan metode limma berhasil mengidentifikasi 733 gen signifikan yang mengalami upregulated dan downregulated. 
- Visualisasi menggunakan volcano plot menunjukkan adanya distribusi gen signifikan yang jelas. 
- Visualisasi heatmap terhadap 50 gen paling signifikan menunjukkan pola ekspresi berbeda secara konsisten antara sampel kanker dan normal. Sampel kanker cenderung mengelompok secara terpisah dari sampel normal, yang menandakan bahwa gen-gen signifikan memiliki potensi sebagai biomarker molekuler dalam membedakan kondisi patologis dan normal. 
- Hasil analisis GO menunjukkan gen yang terekspresi berbeda terlibat dalam proses biologis seperti respons terhadap zat xwnobiotik, metabolisme steroid, respons terhadap hipoksia, regulasi metabolisme molekul kecil, serta proses detoksifikasi. Hal ini mencerminkan adanya perubahan fungsi metabolik dan mekanisme adaptasi sel kanker terhadap lingkungan yang tidak normal. 
- Analisis KEGG menunjukkan gen yang terekspresi berbea berperan dalam berbagai jalur penting, seperti metabolisme obat oleh sitokrom P450, siklus sel, biosintesis hormon steroid, metabolisme retinol, dan karsinogenesis kimia. Jalur tersebut berhubungan dengan proliferasi sel, metabolisme, dan mekanisme perkembangan kanker kolorektal.
