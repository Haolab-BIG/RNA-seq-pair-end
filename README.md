# RNA-seq (pair-end)
For everyone who has a new face to HaoLab, we make this pipeline/protocol to standardize downstream analysis pipeline for RNA-seq (pair-end) data.

## Part I Introduction
### i. Workflow
Here stands an throughout workflow of RNA-seq (pair-end) data analysis.
As illustrated in the figure,
(i) yellow circles represent the steps where commands need to be entered;
(ii) pink dashed rectangular boxes represent the output results after processing at each step.
We will proceed by structuring our workflow according to (ii).
![RNA-seq-pair-end.png](https://github.com/Haolab-BIG/RNA-seq-pair-end/blob/main/RNA-seq-pair-end.png)
### ii. File Structure
Here stands an throughout file structure of RNA-seq (pair-end) data analysis.
* *You can decide what your structure looks like, which makes it more efficient to work.*
```
RNA-seq
├── 1.rawdata
│   ├── iN_CO.rep1_raw_1.fq.gz
│   └── iN_CO.rep1_raw_2.fq.gz
├── 2.FastQC
│   ├── fastqc.log
│   ├── iN_CO.rep1_raw_1_fastqc.html
│   ├── iN_CO.rep1_raw_1_fastqc.zip
│   ├── iN_CO.rep1_raw_2_fastqc.html
│   └── iN_CO.rep1_raw_2_fastqc.zip
├── 3.trim
│   ├── iN_CO.rep1_raw_1.fq.gz_trimming_report.txt
│   ├── iN_CO.rep1_raw_1_val_1_fastqc.html
│   ├── iN_CO.rep1_raw_1_val_1_fastqc.zip
│   ├── iN_CO.rep1_raw_1_val_1.fq.gz
│   ├── iN_CO.rep1_raw_2.fq.gz_trimming_report.txt
│   ├── iN_CO.rep1_raw_2_val_2_fastqc.html
│   ├── iN_CO.rep1_raw_2_val_2_fastqc.zip
│   ├── iN_CO.rep1_raw_2_val_2.fq.gz
│   └── trim-iN_CO.rep1.log
├── 4.StarResult
│   ├── iN_CO.rep1_Aligned.sortedByCoord.out.bam
│   ├── iN_CO.rep1_Aligned.sortedByCoord.out.bam.bai
│   ├── iN_CO.rep1_Log.final.out
│   ├── iN_CO.rep1_Log.out
│   ├── iN_CO.rep1_Log.progress.out
│   ├── iN_CO.rep1_ReadsPerGene.out.tab
│   └── iN_CO.rep1_SJ.out.tab
├── 5.removeDup
│   ├── iN_CO.rep1_marked_dup_metrics.txt
│   ├── iN_CO.rep1_readgroup.bam
│   ├── iN_CO.rep1_removeDup.bam
│   ├── iN_CO.rep1_removeDup.bam.bai
│   ├── read.group_iN_CO.rep1.log
│   └── samflag.log
├── 6.bw
│   ├── bw.log
│   └── iN_CO.rep1_removeDup.bw
├── 7.featureCounts
│   ├── analysis.R
│   ├── analysis.RData
│   ├── diffOutput_IN.log
│   ├── diffOutput_IN.txt
│   ├── featureCounts.log
│   ├── IN.count.txt
│   ├── IN.go.heat.pdf
│   ├── IN.go.pdf
│   ├── IN.heatmap.pdf
│   ├── IN.kegg.heat.pdf
│   ├── IN.kegg.pdf
│   ├── IN.sig.gene.txt
│   ├── PCA_total.pdf
│   ├── reads.summary
│   ├── TotalSample.txt
│   ├── TotalSample.txt.summary
│   ├── total.tpm.txt
│   ├── volcano_bar_IN.pdf
│   └── volcano_IN.pdf
```

### iii. Conda Environment
You can configure a Conda environment named 'RNA-seq-pair-end' using the following code, which includes the essential software for RNA-seq (pair-end) analysis.
```
conda create --name RNA-seq-pair-end python=3.9
conda config --add channels bioconda
conda config --add channels conda-forge
conda install bioconda::trim-galore
conda install bioconda::picard
conda install bioconda::subread
conda install bioconda::deeptools
conda install bioconda::homer
conda install bioconda::star
```

## Part II Generation of Data for Analysis: FASTQ2BAM
In this section, you will convert the raw FASTQ files into BAM files, which can be used for subsequent analysis.
### i. Raw Data Quality Check(QC)
You can perform quality check on the raw data to assess the sequencing quality.
```
mkdir 2.FastQC
fastqc -o 2.FastQC --noextract -f fastq ./1.rawdata/*.fq.gz -t 16 >2.FastQC/fastqc.log 2>&1 &
```
### ii. Trimming After QC
Based on the quality control results, you can perform appropriate trimming on the raw data.
```
mkdir 3.trim
for file in 1.rawdata/*_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_1.fq.gz|\1|')
    echo ${filename}
    trim_galore -q 25 --cores 16 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 \
    --paired ./1.rawdata/${filename}_1.fq.gz ./1.rawdata/${filename}_2.fq.gz \
    -o ./3.trim >3.trim/trim-${filename}.log 2>&1
done
```
### iii. Mapping Clean Data to Genome, Filtering and Remove Duplicates
Align the cleaned data obtained in the previous step to the reference genome of the corresponding species. After filtering for high-quality sequences and removing duplicates, the resulting data will be prepared for further analysis.

#### 1. Mapping, Filtering
Align the quality-controlled FASTQ files to the reference genome, filter for high-quality aligned sequences, and perform sorting.

```
mkdir 4.StarResult
IndexPath=
TmpPath=
for file in 1.rawdata/*_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_1.fq.gz|\1|')
    echo ${filename}
    STAR --runMode alignReads --genomeDir ${IndexPath} \
    --outFileNamePrefix ./4.StarResult/${filename}_ \
    --outSAMattributes All --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ./3.trim/${filename}_1_val_1.fq.gz ./3.trim/${filename}_2_val_2.fq.gz \
    --runThreadN 20 --readFilesCommand zcat --outTmpDir ${TmpPath}/${filename} \
    --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 \
    –-genomeLoad LoadAndKeep >4.StarResult/Star_${filename}.log 2>&1
done
STAR --genomeLoad Remove --genomeDir ${indexPath}
```
#### 2. Add read group
Add read group for each sample.
```
mkdir 5.removeDup
for file in 1.rawdata/*_raw_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    echo ${filename}_Aligned.sortedByCoord.out.bam
    samtools index -@ 16 4.StarResult/${filename}_Aligned.sortedByCoord.out.bam
    picard AddOrReplaceReadGroups \
    I=4.StarResult/${filename}_Aligned.sortedByCoord.out.bam \
    O=5.removeDup/${filename}_readgroup.bam \
    RGID=${filename} \
    RGLB=${filename} \
    RGPL=illumina \
    RGPU=${filename} \
    RGSM=${filename} > 5.removeDup/read.group_${filename}.log 2>&1 &
done
```

#### 3. Remove Duplicates
Remove duplicate sequences.
```
for file in 1.rawdata/*_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    picard MarkDuplicates I=5.removeDup/${filename}_readgroup.bam O=5.removeDup/${filename}_removeDup.bam M=5.removeDup/${filename}_marked_dup_metrics.txt REMOVE_DUPLICATES=true READ_NAME_REGEX=null >5.removeDup/MarkDuplicates_${filename}.log 2>&1 &
done
for file in 1.rawdata/*_raw_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    echo ${filename}_removeDup.bam >>5.removeDup/samflag.log
    samtools flagstat -@ 10 5.removeDup/${filename}_removeDup.bam >> 5.removeDup/samflag.log
done
```
## Part III Visualization on Genome Browser: BAM2BW
Convert the alignment results into files with lower resolution but smaller size for easier browsing across the genome; simultaneously, assess the correlation among the experimental samples.

### i. Visualization on Genome Browser
You can visualize these results in the UCSC genome browser or IGV genome browser.
```
mkdir 6.bw
for file in 1.rawdata/*_raw_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    echo ${filename}_removeDup.bam >> 7.bw/bw.log
    samtools index -@ 16 5.removeDup/${filename}_removeDup.bam
    bamCoverage -b 5.removeDup/${filename}_removeDup.bam -o 7.bw/${filename}_removeDup.bw  -p 30 --binSize 1 --scaleFactor 1 --normalizeUsing CPM >>  7.bw/bw.log 2>&1
done
```
## Part IV Analysis of Gene Expression
In this section, you will acquire gene expression data and perform a series of analyses based on it.

### i. Quantification for genomic features: BAM2TAB
Count mapped reads for genomic features such as genes.
```
mkdir 7.featureCounts
featureCounts -a ./Index/gencode.v47.annotation.gtf -s 2 -p --countReadPairs -B -T 30  --ignoreDup -o ./7.featureCounts/TotalSample.txt 5.removeDup/*_removeDup.bam  > ./7.featureCounts/featureCounts.log 2>&1 &
```
### ii.Identify differentially expressed genes

#### 1. Scale the gene expression
Calculate TPM(TPM, Transcripts Per Kilobase of exon model per Million mapped reads) through gene length and sequencing depth.
```
featureCountsPath <-
setwd(featureCountsPath)
total.counts <- read.table("TotalSample.txt", sep = "\t", header= T, quote = "")
total.counts <- total.counts[,-c(2,3,4,5)]
total.counts$Geneid <- sub("\\..*$","",total.counts$Geneid)
colnames(total.counts) <- gsub("X5.removeDup.|_removeDup.bam", "", colnames(total.counts))
total.feature <- total.counts
rownames(total.feature) <- total.feature$Geneid
total.feature <- total.feature[,-1]
total.feature$kb <- total.feature$Length/1000
total.rpk <- total.feature[,c(2:(ncol(total.feature)-1))] / total.feature$kb
total.tpm <- as.data.frame(t(t(total.rpk)/colSums(total.rpk) * 1000000)
```
#### 2. Observe batch effect
Batch effects between samples can be identified, allowing for either correction or selective sample screening as needed.
```
pca.info <- fast.prcomp(total.tpm)
pca.data <- data.frame(sample = rownames(pca.info$rotation),pca.info$rotation)
pc_contribution <- (pca.info$sdev^2) / sum(pca.info$sdev^2)
pc1_contribution <- round(pc_contribution[1] * 100, 2)
pc2_contribution <- round(pc_contribution[2] * 100, 2)
p <- ggscatter(pca.data, x = "PC1", y = "PC2", color = "sample") +
  scale_color_manual(values = c(
    "231_CO"="#8c510a","231_CO.rep1"="#bf812d","231_CO.rep2"="#dfc27d","231_CO.rep3"="#f6e8c3",
    "231_MO"="#01665e","231_MO.rep1"="#35978f","231_MO.rep2"="#80cdc1","231_MO.rep3"="#c7eae5",
    "iN_CO"="#762a83","iN_CO.rep1"="#9970ab","iN_CO.rep2"="#c2a5cf","iN_CO.rep3"="#e7d4e8",
    "iN_MO"="#2166ac","iN_MO.rep1"="#4393c3","iN_MO.rep2"="#92c5de","iN_MO.rep3"="#d1e5f0"
  )) +
  theme_base() +
  labs(
    x = paste("PC1 (", pc1_contribution, "%)", sep=""),
    y = paste("PC2 (", pc2_contribution, "%)", sep="")
  )
ggsave("PCA_total.pdf", plot = p, device = "pdf", width = 8, height = 4, path = featureCountsPath)
IN.count <- total.counts[,c(1,12:14,16:18)]  #Geneid and selected samples
write.table(IN.count,"IN.count.txt",sep = "\t",col.names = T,row.names = F,quote = F)
```
#### 3. Calculate the difference between two groups
Calculate logFC and p value for each gene between two groups using read counts.

```
getDiffExpression.pl 7.featureCounts/IN.count.txt IN_CO IN_CO IN_CO IN_MO IN_MO IN_MO -edgeR > 7.featureCounts/diffOutput_IN.txt 2> 7.featureCounts/diffOutput_IN.log &
```
#### 4. Filter out differentially expressed genes
Set filtering criterion, and obtain differentially expressed genes.
```
IN.homer <-  read.table("diffOutput_IN.txt",sep = "\t", header = T, quote = "")
colnames(IN.homer)[1:7] <- colnames(IN.count)
colnames(IN.homer)[8:10] <- c("logFC","p","p.adj")
IN.up <- IN.homer[IN.homer$p<0.05 & IN.homer$logFC>1,]
IN.down <- IN.homer[IN.homer$p<0.05 & IN.homer$logFC< -1,]
IN.notsig <- IN.homer[!IN.homer$Geneid%in%IN.up$Geneid & !IN.homer$Geneid%in%IN.down$Geneid,]
```
Draw volcano and heat maps
```
p <-ggplot() +
  geom_point(data = IN.notsig, aes(x = logFC, y = -log10(p)), size = 1, col = "grey") +
  geom_point(data = IN.up, aes(x = logFC, y = -log10(p)), col = "#c51b7d", size = 1) +
  geom_point(data = IN.down, aes(x = logFC, y = -log10(p)), col = "#4d9221", size = 1) +
  theme_bw() +
  ggtitle("IN_CO VS IN_MO")+
  ylim(0, 25) + 
  xlim(-7.5,7.5) +
  ylab(expression("-log"[10]*" (p)")) +
  xlab(expression("log"[2]*"FC")) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggsave("volcano_IN.pdf", plot = p, device = "pdf", width = 3, height = 3, path = featureCountsPath)
IN.table <- data.frame(Group = c("up","down") ,number = c(nrow(IN.up),nrow(IN.down)))
p <- ggplot(IN.table, aes(x = Group, y = number, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = number, y = number / 2), hjust = 1.5, color="white") +
  scale_fill_manual(values = c("up" = "#c51b7d", "down" = "#4d9221")) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = NULL, x = NULL, y = NULL)
ggsave("volcano_bar_IN.pdf", plot = p, device = "pdf", width = 1, height = 0.5, path = featureCountsPath)
IN.sig.tpm <- total.tpm[rownames(total.tpm) %in% c(IN.up$Geneid, IN.down$Geneid),
                          intersect(colnames(total.tpm), colnames(IN.count))]
IN.sig.tpm.nor <- as.data.frame(t(scale(t(IN.sig.tpm))))
p<-pheatmap(IN.sig.tpm.nor, 
            cluster_cols = F,
            cluster_rows = T,
            clustering_distance_rows = "euclidean", 
            clustering_method = "ward.D2",
            cutree_rows = 2,
            #annotation_row=diff.sig.clusters,
            #annotation_colors = annotation_colors,
            show_rownames = F,
            show_colnames = TRUE)
ggsave("IN.heatmap.pdf", plot = p, device = "pdf", width = 4, height = 8, path = featureCountsPath)
```
#### 5.Functional analysis
GO (GO, Gene Ontology) Enrichment Analysis
```
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
IN.sig.gene <- rbind(IN.up,IN.down)
IN.sig.gene <- IN.sig.gene %>% select(Geneid, logFC, p, p.adj)
row.names(IN.sig.gene) <- IN.sig.gene$Geneid
IN.sig.gene$symbol = mapIds(x= org.Hs.eg.db, 
                            keys = IN.sig.gene$Geneid,
                            keytype="ENSEMBL",
                            column ="SYMBOL",
                            multiVals = "first")
IN.sig.gene$entriz = mapIds(x= org.Hs.eg.db, 
                            keys = IN.sig.gene$Geneid,
                            keytype="ENSEMBL",
                            column ="ENTREZID",
                            multiVals = "first")
write.table(IN.sig.gene,"IN.sig.gene.txt",sep="\t",col.names = T,row.names = F,quote = F)
```
```
IN.sig.gene.old <- na.omit(IN.sig.gene)
logFC <- IN.sig.gene.old$logFC
names(logFC) <- IN.sig.gene.old$symbol
library(clusterProfiler)
GOenrich = enrichGO(gene = IN.sig.gene.old$entriz, 
                      OrgDb = org.Hs.eg.db,
                      keyType ="ENTREZID",
                      pAdjustMethod ='BH',
                      ont = "ALL",
                      pvalueCutoff=0.05,
                      qvalueCutoff=0.05,
                      readable =T)
top10_results <- GOenrich
GOenrich.top10 <- as.data.frame(GOenrich@result) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)
top10_results@result <- GOenrich.top10
pdf(file = "IN.go.heat.pdf", width =12, height = 4)
heatplot(top10_results, foldChange=logFC)
dev.off()
GOenrich.top10$logp.adj <- -log10(GOenrich.top10$p.adjust)
GOenrich.plot <- ggplot(data=GOenrich.top10, aes(x=Description, y=logp.adj )) +
  geom_bar(stat = "identity", fill = "#d33682", width = 0.8) + 
  coord_flip() + 
  theme_test() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlab("GO term") + 
  theme(axis.text=element_text(color="black")) +
  labs(title = "TOP 10 enriched GO term", y = expression("-log"[10]*" (p.adj)"),  x = "")
ggsave("IN.go.pdf", plot = GOenrich.plot, width = 5, height = 5)

```
KEGG (KEGG, Kyoto Encyclopedia of Genes and Genomes) Enrichment Analysis
```
KEGGenrich = enrichKEGG(gene =IN.sig.gene.old$entriz, 
                          keyType = "kegg",
                          pAdjustMethod ='BH',
                          organism= "hsa",
                          qvalueCutoff =0.05,
                          pvalueCutoff=0.05)
top10_results <- KEGGenrich
KEGGenrich.top10 <- KEGGenrich[order(KEGGenrich$p.adjust,decreasing = F)[1:10],]
top10_results@result <- KEGGenrich.top10
top10_results <- setReadable(top10_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
pdf(file = "IN.kegg.heat.pdf", width =12, height = 4)
heatplot(top10_results, foldChange=logFC)
dev.off()
KEGGenrich.top10$logp.adj <- -log10(KEGGenrich.top10$p.adjust)
KEGGenrich.plot <- ggplot(data=KEGGenrich.top10, aes(x=Description, y=logp.adj )) +
  geom_bar(stat = "identity", fill = "#d33682", width = 0.8) + 
  coord_flip() + 
  theme_test() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlab("KEGG pathway") + 
  theme(axis.text=element_text(color="black")) +
  labs(title = "TOP 10 enriched KEGG pathway", y = expression("-log"[10]*" (p.adj)"),  x = "")
ggsave("IN.KEGG.pdf", plot = KEGGenrich.plot, width = 8, height = 5)
```






