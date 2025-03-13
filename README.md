# RNA-seq (pair-end)
For everyone who has a new face to HaoLab, we make this pipeline/protocol to standardize downstream analysis pipeline for RNA-seq (pair-end) data.

## Part I Introduction
### i. Workflow
Here stands an throughout workflow of RNA-seq (pair-end) data analysis.
As illustrated in the figure,
(i) yellow circles represent the steps where commands need to be entered;
(ii) pink dashed rectangular boxes represent the output results after processing at each step.
We will proceed by structuring our workflow according to (ii).

### ii. File Structure
Here stands an throughout file structure of RNA-seq (pair-end) data analysis.
* *You can decide what your structure looks like, which makes it more efficient to work.*

### iii. Conda Environment
You can configure a Conda environment named 'RNA-seq-pair-end' using the following code, which includes the essential software for RNA-seq (pair-end) analysis.
```
conda create --name RNA-seq-pair-end python=3.9
conda config --add channels bioconda conda-forge

conda install trim-galore picard subread deeptools homer


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
for file in 1.rawdata/*_raw_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    echo ${filename}
    trim_galore -q 25 --cores 16 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 \
    --paired ./1.rawdata/${filename}_raw_1.fq.gz ./1.rawdata/${filename}_raw_2.fq.gz \
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
for file in 1.rawdata/*_raw_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|1.rawdata/(.*)_raw_1.fq.gz|\1|')
    echo ${filename}
    STAR --runMode alignReads --genomeDir ${indexPath} \
    --outFileNamePrefix ./4.StarResult/${filename}_ \
    --outSAMattributes All --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ./3.trim/${filename}_raw_1_val_1.fq.gz ./3.trim/${filename}_raw_2_val_2.fq.gz \
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
for file in 1.rawdata/*_raw_1.fq.gz; do
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
mkdir 7.bw
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
mkdir 6.featureCounts
featureCounts -a ./Index/gencode.v47.annotation.gtf -s 2 -p --countReadPairs -B -T 30  --ignoreDup -o ./6.featureCounts/TotalSample.txt 5.removeDup/*_removeDup.bam  > ./6.featureCounts/featureCounts.log 2>&1 &
```
### ii.Identify differentially expressed genes

#### 1. Observe batch effect
Batch effects between samples can be identified, allowing for either correction or selective sample screening as needed.
```

```
#### 2. Scale the gene expression
Calculate TPM(TPM, Transcripts Per Kilobase of exon model per Million mapped reads) through gene length and sequencing depth.

#### 3. Calculate the difference between two groups
Calculate logFC and p value for each gene between two groups using read counts.
```
getDiffExpression.pl 6.featureCounts/IN.count.txt IN_CO IN_CO IN_CO IN_MO IN_MO IN_MO -edgeR > 6.featureCounts/diffOutput_IN.txt 2> 6.featureCounts/diffOutput_IN.log &
```
#### 4. Filter out differentially expressed genes
Set filtering criterion, and obtain differentially expressed genes.


Draw volcano and heat maps


#### 5.Functional analysis
GO (GO, Gene Ontology) Enrichment Analysis


KEGG (KEGG, Kyoto Encyclopedia of Genes and Genomes) Enrichment Analysis










