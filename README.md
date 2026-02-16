RNA-Seq Differential Expression Analysis Pipeline

This repository contains a complete RNA-seq analysis workflow from raw FASTQ files to differential gene expression, visualization, and interpretation across multiple stress conditions.

Project Overview
This pipeline performs:

- Quality control of raw reads
- Adapter trimming
- Genome alignment
- Gene quantification
- Differential expression analysis
- Visualization (PCA, Volcano, Heatmap, MA plots)

Organism: *Cyamopsis tetragonoloba*  
Conditions analyzed:
- Control
- Heat
- Drought
- Salinity

Each condition contains **3 biological replicates**

RNA-Seq Analysis Pipeline

1. Download FASTQ files
```bash
fasterq-dump --split-files SRRXXXXXXX
```

2. Quality Control
Falco: Galaxy

3. Trimming Adaptors: Trimmomatic
```bash
Download trimmomatic: github
java -jar /mnt/d/gaurgram/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 6 \
control_R1_rep1.fastq  control_R2_rep1.fastq \
control_R1_rep1_paired.fastq control_R1_rep1_unpaired.fastq \
control_R2_rep1_paired.fastq control_R2_rep1_unpaired.fastq \
ILLUMINACLIP:/mnt/d/gaurgram/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Input: forward and reverse R1/R2 files
Output: paired and unpaired R1/R2 files
```

4. Index Reference Genome
```bash
hisat2-build ref_genome.fa index_prefix           [No gtf file is avaliable for Cyamopsis tetragonoloba] 
```

5. Alignment and sorting
```bash 
hisat2 -p 4 --dta -x index_prefix \
-1 control_R1_rep2.fastq -2 control_R2_rep2.fastq \
| samtools sort -o control_rep2.sorted.bam

Input: R1/R2 fastq files
Output: sorted.bam files
```

6. Quantify gene expression
```bash
featureCounts -T 4 \
-p -B -C \
-t exon \
-g gene_id \
-a merged.gtf \
-o feature_counts.txt \
*.sorted.bam
```

7. Load Counts matrix in R
```r
counts <- read.delim("feature_counts.txt",
                     comment.char="#",
                     row.names=1)
```

8. Metadata
```r
colData <- data.frame(
 condition=c(
 "control","control","control",
 "heat","heat","heat",
 "drought","drought","drought",
 "salinity","salinity","salinity"
))
rownames(colData) <- colnames(counts)
```
<img width="321" height="246" alt="image" src="https://github.com/user-attachments/assets/9ba864c5-fb03-49e0-b74f-35e4958e08ae" />

9. DESeq2 Analysis
```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
 countData=counts,
 colData=colData,
 design=~condition)

dds <- dds[rowSums(counts(dds))>10,]
dds <- DESeq(dds)
```

10. Differential Expression
```r
res_heat <- results(dds, contrast=c("condition","heat","control"))
res_drought <- results(dds, contrast=c("condition","drought","control"))
res_salinity <- results(dds, contrast=c("condition","salinity","control"))
```

11. Normalization
```r
vsd <- vst(dds, blind=FALSE)
vst_mat <- assay(vsd)
```
12. Visualization
    
PCA Plot
```r
pcaData <- plotPCA(vsd,intgroup="condition",returnData=TRUE)

ggplot(pcaData,aes(PC1,PC2,color=condition))+
geom_point(size=4)+
theme_bw()
```
<img width="613" height="351" alt="pca_ggplot2" src="https://github.com/user-attachments/assets/13e7a249-d886-4380-996f-e09229546a03" />

Volcano Plot
```r
res_df <- as.data.frame(res_drought)
res_df$gene <- rownames(res_df)

res_df$diffexpressed <- "Not Significant"
res_df$diffexpressed[res_df$log2FoldChange>=1 & res_df$padj<0.05] <- "Up"
res_df$diffexpressed[res_df$log2FoldChange<=-1 & res_df$padj<0.05] <- "Down"

ggplot(res_df,aes(log2FoldChange,-log10(padj),color=diffexpressed))+
geom_point()
```
<img width="4800" height="1800" alt="combined_volcano_stress_conditions" src="https://github.com/user-attachments/assets/46d3549a-505d-453b-bf77-1ac0c1ebbfb4" />

Heatmap
```r
shared_genes <- Reduce(intersect,
 list(genes_drought,genes_heat,genes_salinity))

pheatmap(heatmap_mat)
```
![shared_DEG_heatmap](https://github.com/user-attachments/assets/f86fe932-7b3e-449b-bf27-607b011a9551)

MA Plots
```r
plotMA(res_drought)
plotMA(res_heat)
plotMA(res_salinity)
```
![MA_plots_all_stresses](https://github.com/user-attachments/assets/76b2c544-7541-4f59-8694-5ac9ce3a6787)




