# This script is generated to analyse differential gene expression between human peritoneal dialysis effluent cells treated with or without alanyl-glutamine supplemented peritoneal dialysis fluid in a randomized controlled cross-over pilot study
#Source Array express (http://www.ebi.ac.uk/arrayexpress/) 
#ArrayExpressAccession:	E-MTAB-5462
#Illumina single end (1x50 bp) sequencing. Each read is 50bp.
#This script was produced following the RsubreadUserguide (https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/Rsubread.pdf)
## Trimmed Truseq adapter Index1/18 (GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTAT) from P1-AG.fastq.gz file and retained sequences with Phred quality >= 33.
#Built human reference genome index file and aligned short read sequences in fastq file to human hg19 reference genome using subread aligner
#Generated .bam files using Rsubread align
#Assigned mapped sequencing reads to hg19 Entrez genes
#generated featureCounts
#Generated 'countmatrix' (Countdata) using DESeq2 package.
#estimating dispersions

##WARNING: 
#Warning message:Accurate differential expression analysis requires replication AND without Experimental Replicates we can not calculae the the biological coeffiecient of variation as variation requires sample mean. Hence this exercize is only for data exloration. The differences between samples can be viewed in terms of fold change but should not rely on the strength of the significance. 

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
#load the workspace
setwd("/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR")
load(file="PR_RNAseq-Test_Oxford.RData")
library(DESeq2)
counts<-feature_counts$counts
head(counts)
dim(counts)

# Convert to matrix
countdata <- as.matrix(counts)
head(countdata)

#Define design data frame based on the experimental design.

Treat<-factor(c("AG", "Non-AG"), levels=c("Non-AG", "AG"))
design<-data.frame(Treat)


# Create a coldata frame using design data frame and create the dds (DESeqDataSet) using DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), Treat))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~Treat)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Differential expression 
de <- results(dds)
table(de$padj<0.05)

## Order by adjusted p-value
de <- de[order(de$padj), ]

## Merge with normalized count data
dedata <- merge(as.data.frame(de), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(dedata)[1] <- "Gene"
head(dedata)
## Write differentially expressed genes list to table
write.csv(dedata, file="PR_RNAseq-Test_Oxford_DGE_P1-AG_vs_P1.csv")


