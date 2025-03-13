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

## Part II Generation of Data for Analysis: FASTQ2BAM
In this section, you will convert the raw FASTQ files into BAM files, which can be used for subsequent analysis.
### i. Raw Data Quality Check(QC)
You can perform quality check on the raw data to assess the sequencing quality.

### ii. Trimming After QC
Based on the quality control results, you can perform appropriate trimming on the raw data.

### iii. Mapping Clean Data to Genome, Filtering and Remove Duplicates
Align the cleaned data obtained in the previous step to the reference genome of the corresponding species. After filtering for high-quality sequences and removing duplicates, the resulting data will be prepared for further analysis.

#### 1. Mapping, Filtering
Align the quality-controlled FASTQ files to the reference genome, filter for high-quality aligned sequences, and perform sorting.

#### 2. BAM Indexing
Index the sorted files.

#### 3. Add read group
Add read group for each sample.

#### 4. Remove Duplicates
Remove duplicate sequences.

## Part III Visualization on Genome Browser: BAM2BW
Convert the alignment results into files with lower resolution but smaller size for easier browsing across the genome; simultaneously, assess the correlation among the experimental samples.

### i. Visualization on Genome Browser
You can visualize these results in the UCSC genome browser or IGV genome browser.


## Part IV 
Count mapped reads for genomic features such as genes.










