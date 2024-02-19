#Tutorial_Module3_Part4_edgeR.R

#Malachi Griffith, mgriffit[AT]genome.wustl.edu
#Obi Griffith, ogriffit[AT]genome.wustl.edu
#The Genome Institute, Washington Univerisity School of Medicine
#R tutorial for CBW - Informatics for RNA-sequence Analysis

#######################
# Loading Data into R #
#######################
library(edgeR)
#Set working directory where output will go
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/counts/edgeR/transcript_level/SaltVsNosalt"
setwd(working_dir)

#Read in gene mapping
mapping=read.csv("transcriptinfo.csv", header=FALSE, stringsAsFactors=FALSE, row.names=1)

# Read in count matrix
dat=read.csv("SaltVsNosalt_transcript_count_matrix.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

#The last 5 rows are summary data, remove
rawdata=dat

# Set column names (optional, but this helps keep the conditions straight)
#  colnames(rawdata) <- c("Bac_rep1","Bac_rep2","Bac_rep3","Metabact_rep1","Metabact_rep2","Metabact_rep3")
#colnames(rawdata) <- c("Salt_rep1","Salt_rep2","Salt_rep3","Metabact_rep1","Metabact_rep2","Metabact_rep3")
 # colnames(rawdata) <- c("Salt_rep1","Salt_rep2","Salt_rep3","Bac_rep1","Bac_rep2","Bac_rep3")
colnames(rawdata) <- c("Salt_rep1","Salt_rep2","Salt_rep3","NoSalt_rep1","NoSalt_rep2","NoSalt_rep3")

# Check dimensions
dim(rawdata)

# Require at least 25% of samples to have count > 25
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)

#################
# Running edgeR #
#################
library("limma", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
# load edgeR
library('edgeR')

# make class labels
class <- factor( c( rep("Salt",3), rep("Nosalt",3) ))

# Get common gene names
genes=rownames(rawdata)
gene_names=mapping[genes,1]


# Make DGEList object
y <- DGEList(counts=rawdata, genes=genes, group=class)
nrow(y)

# TMM Normalization
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Differential expression test
et <- exactTest(y)

# Print top genes
topTags(et)

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
head(et$table)
# Output DE genes
# Matrix of significantly DE genes
mat <- cbind(
  genes,gene_names, et$table$PValue,
  sprintf('%0.3f',log10(et$table$PValue)),
  sprintf('%0.3f',et$table$logFC)
 )[as.logical(de),]
head(mat)

colnames(mat) <- c("Gene", "Gene_Name", "PValue", "Log10_Pvalue", "Log_fold_change")

# Order by log fold change
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]
dim(mat)
# Save table
write.table(mat, file="DE_transcripts_SaltVsNosalt.txt", quote=FALSE, row.names=FALSE, sep="\t")

#To exit R type the following
#quit(save="no")
