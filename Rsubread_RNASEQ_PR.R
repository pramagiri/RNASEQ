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
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("Rsubread")
biocLite("edgeR")
biocLite("limma")
biocLite("Homo.sapiens")
install.packages("xlsx")
install.packages("locfit")

library(Rsubread)
library(limma)
library(edgeR)
setwd("/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR")
buildindex(basename="hg19",reference="hg19.fa")

fasta_files1<-list.files(path="/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR/FQ", pattern="trimmed_P1.fastq.gz")
align(index="hg19",readfile1=fasta_files1, output_file="alignP1.bam")
fasta_files2<-list.files(path="/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR/FQ", pattern="trimmed_P1-AG.fastq.gz")
align(index="hg19",readfile1=fasta_files2, output_file="alignP1-AG.bam")

#Assigned mapped sequencing reads to hg19 Entrez genes
bam_files<-list.files(path="/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR/FQ/", pattern=".bam$")
bam_files<-paste("/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR/",bam_files, sep="")

annotation<-read.table("/media/pradeep/NGSPR_HOME/RNASeq-Test_Oxford_PR/hg19_featureCount.txt", header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!" )
saf<-annotation
feature_counts<-featureCounts(files=bam_files, annot.ext=saf, allowMultiOverlap= FALSE, isPairedEnd=FALSE)

#Replace long sample names with short ones
head(feature_counts$counts)
colnames(feature_counts$counts)<-sub("X.media.pradeep.NGSPR_HOME.RNASeq.Test_Oxford_PR.FQ.trimmed", "", colnames(feature_counts$counts))
colnames(feature_counts$counts)<-sub(".bam", "", colnames(feature_counts$counts))
colnames(feature_counts$counts)<-sub("_", "", colnames(feature_counts$counts))
head(feature_counts$counts)
save.image(file="PR_RNAseq-Test_Oxford.RData")
