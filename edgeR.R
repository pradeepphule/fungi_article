library(supraHex)

# Load or install packages (i.e. edgeR) specifically used in this demo
for(pkg in c("edgeR")){
  if(!require(pkg, character.only=T)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg)
    lapply(pkg, library, character.only=T)
  }
}
# Import data
## import RNA-seq counts data file RNAseq_counts.txt
RNAseq_counts <- read.csv(file = "/home/pradeep/RNA_HOME/gene_count_matrix.csv", header = TRUE, row.names = 1)

## import RNA-seq counts data file RNAseq_geneinfo.txt
RNAseq_geneinfo <- read.csv(file = "/home/pradeep/RNA_HOME/RNAseq_geneinfo.csv", header = TRUE, row.names = 1) 


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Identify differentially expressed genes using the package 'edgeR'
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# (I) Create a DGEList object (edgeR's container for RNA-seq count data)
d_obj <- DGEList(counts=RNAseq_counts, genes=RNAseq_geneinfo)

## In edgeR, it is recommended to remove genes/features without at least 1 read/count per million (known as 'cpm') in n of the samples, where n is the size of the smallest group of replicates. In this case, n=3 but we set it to be 6 as the strictest filtering
cpms <- edgeR::cpm(d_obj$counts)
keep <- rowSums(cpms>=1) >= 6
d <- d_obj[keep,]

## Reset the library sizes
d$samples$lib.size <- colSums(d$counts)

## Estimate normalization factors using:
d <- calcNormFactors(d)

## Define the design matrix
cnames <- colnames(d)
cellline <- gsub("_.*", "", cnames, perl=T)
## treatment <- gsub(".*_", "", cnames, perl=T)
targets <- data.frame(sample=cnames, cellline=cellline)
design <- model.matrix(~cellline, targets)
design

## Inspect the relationships between samples using a multidimensional scaling (MDS) plot to show the relationship between all pairs of samples
plotMDS(d, labels=colnames(d), col=c("red","darkgreen","blue","yellow")[factor(targets$cellline)], xlim=c(-5,5), ylim=c(-5,5))
### Note: this inspection clearly shows the variances between cell lines, calling for paired design test.

# (II) Do dispersion estimation and GLM (Generalized Linear Model) fitting
# Estimate the overall dispersion to get an idea of the overall level of biological variablility


## Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood:
 d2 <- estimateGLMTrendedDisp(d, design)
d2 <- estimateGLMTagwiseDisp(d2, design)

## Given the design matrix and dispersion estimates, fit a GLM to each feature:
f <- glmFit(d2, design)

# (III) Perform a likelihood ratio test
## Specify the difference of interest: DEX vs CON
contrasts <- rbind(c(0,0,0,1))

## Prepare the results
logFC <- matrix(nrow=nrow(d2), ncol=1)
PValue <- matrix(nrow=nrow(d2), ncol=1)
FDR <- matrix(nrow=nrow(d2), ncol=1)
tmp <- c("DEX_CON")
colnames(logFC) <- paste(tmp, '_logFC', sep='')
colnames(PValue) <- paste(tmp, '_PValue', sep='')
colnames(FDR) <- paste(tmp, '_FDR', sep='')
rownames(logFC) <- rownames(PValue) <- rownames(FDR) <- rownames(d2)

## Perform the test, calculating P-values and FDR (false discovery rate)
for(i in 1:nrow(contrasts)){   
  lrt <- glmLRT(f, contrast=contrasts[i,])
  tt <- topTags(lrt, n=nrow(d2), adjust.method="BH", sort.by="none")
  logFC[,i] <- tt$table$logFC
  PValue[,i] <- tt$table$PValue
  FDR[,i] <- tt$table$FDR
}

## MA plots for RNA-seq data
lrt <- glmLRT(f, contrast=contrasts[1,])
#lrt <- glmLRT(f, coef=4)
summary(de <- decideTestsDGE(lrt, adjust.method="BH", p.value=0.05))

detags <- rownames(d2)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
### As you have seen, it plots the log-fold change (i.e., the log ratio of normalized expression levels between two experimental conditions (i.e. DEX vs CON) against the log counts per million (CPM). Those genes selected as differentially expressed (with a 5% false discovery rate) are highlighted as red dots

# (IV) Output edgeR results
## log counts per million (CPM) for each sample
cpms <- edgeR::cpm(d2,log=T,prior.count=2)
colnames(cpms) <- paste(colnames(d2), '_CPM', sep='')

## log ratio between the treatment vs the control (for each cell line)
odd_indexes <- seq(1,ncol(cpms),2)
even_indexes <- seq(2,ncol(cpms),2)
logFC_cpm <- cpms[,even_indexes] - cpms[,odd_indexes]
colnames(logFC_cpm) <- gsub("_CPM", "_CPM_logFC", colnames(logFC_cpm), perl=T)

## write into the file 'RNAseq_edgeR.txt'
out <- data.frame(EnsemblGeneID=rownames(d2$genes), d2$genes, cpms, logFC_cpm, logFC, PValue, FDR)
write.table(out, file="RNAseq_edgeR.txt", col.names=T, row.names=F, sep="\t", quote=F)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Analyse differentially expressed genes using the package 'supraHex'
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

# (I) Load the package and select differentially expressed genes identified by edgeR data
library(supraHex)
select_flag <- FDR<0.05
data <- logFC_cpm[select_flag, ]
## check the data dimension: how many genes are called significant
dim(data)
# (II) Train the supra-hexagonal map with input data only
sMap <- sPipeline(data, xdim=8, algorithm="sequential")
visHexMulComp(sMap, title.rotate=5, colormap="darkgreen-lightgreen-lightpink-darkred")
