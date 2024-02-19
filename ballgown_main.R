library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data from a file we saved in the current working directory
#/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_Bac_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_Metabact_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Bac_vs_Metabact_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/ref_guided/ballgown/Salt_vs_NoSalt_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/ref_guided/ballgown/Salt_vs_Bac_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/ref_guided/ballgown/Salt_vs_Metabact_path.csv
#/home/pradeep/RNA_HOME/fusarium/stringtie/ref_guided/ballgown/Bac_vs_Metabact_path.csv


pheno_data = read.csv("/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_path.csv")

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path),pData=pheno_data, meas='all')

# Display a description of this object
bg

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
bg_filt
# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
dim(results_transcripts)
dim(results_genes)
head(results_transcripts)
head(results_genes)
# Load all attributes including gene name
bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

## Add gene name
##results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
##results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt),
##                               geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
##results_transcripts <- data.frame(transcriptNames=ballgown::transcriptNames(bg_filt),
##                               transcriptIDs=ballgown::transcriptIDs(bg_filt), results_transcripts)
##  results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt),
##               geneIDs=ballgown::geneIDs(bg_filt),results_transcripts)
## results_transcripts = merge(results_transcripts,bg_table,by.x=c("id"),by.y=c("gene_id"))


results_transcripts <- data.frame(results_transcripts,log_fc = log(results_transcripts
                                                          $fc),Log10_Pvalue = log10(results_transcripts$pval))
results_genes <- data.frame(results_genes,log_fc = log(results_genes
                  $fc),Log10_Pvalue = log10(results_genes$pval))
dim(results_transcripts)
dim(results_genes)

results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
results_transcripts = merge(results_transcripts,bg_table,by.x=c("id"),by.y=c("t_id"))
#results_transcripts = merge(results_transcripts,bg_table$t_name,by.x=c("id"),by.y=c("t_id"))
dim(results_transcripts)
dim(results_genes)
head(results_transcripts)
head(results_genes)
## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)


# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts,results_transcripts$pval<0.05)
sig_genes = subset(results_genes,results_genes$pval<0.05)
dim(sig_transcripts)
dim(sig_genes)

head(sig_transcripts)
head(sig_genes)
## get the p-values and calculate the scores thereupon
sig_genes_pval <- sig_genes$pval
## look at the distribution of p-values
hist(sig_genes_pval)
## get the p-values and calculate the scores thereupon
sig_transcripts_pval <- sig_transcripts$pval
## look at the distribution of p-values
hist(sig_transcripts_pval)
## Sort results from log fold change
sig_transcripts <-  arrange(sig_transcripts, desc(log_fc))
sig_genes <-  arrange(sig_genes, desc(log_fc))


# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_transcript_results_sig.tsv",sep="\t")
write.table(sig_genes,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_gene_results_sig.tsv",sep="\t")
# Output the signifant gene results to a pair of tab delimited files
write.table(results_transcripts,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_transcript_results.tsv",sep="\t")
write.table(results_genes,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_gene_results.tsv",sep="\t")

library(limma)
limmaUsersGuide()
?topSplice()
#fit a linear model from exon level read count data
fit <- lmFit(eexpr(bg_filt), c(0,0,0,1,1,1))
dim(fit)
transcript_id_by_exon = indexes(bg_filt)$e2t$t_id[match(unique(indexes(bg_filt)$e2t$e_id), indexes(bg_filt)$e2t$e_id)]
transcript_id_by_gene = indexes(bg_filt)$t2g$t_id #transcript/gene mapping
geneID = indexes(bg_filt)$t2g$g_id[match(transcript_id_by_exon, transcript_id_by_gene)]
# collect results from diffsplice
ex <- diffSplice(fit, geneID)
write.table(ex,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_diffSplice_counts.tsv",sep="\t")

##  The F-statistic tests for any differences in exon usage between experimental conditions.
topsplice_simes<-topSplice(ex, test="simes", n=20,FDR=1)
write.table(topsplice_simes,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_topsplice_simes_counts.tsv",sep="\t")
## i will prefer f test results as it takes into count exon numbers as well 
##if exon no. are more from same gene then there is more possibility of differentially spliceing from same gene in differnet conditions.
topsplice_ftest<-topSplice(ex, test="F", n=20,FDR=1)
write.table(topsplice_ftest,"/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/Salt_vs_NoSalt_topsplice_ftest_counts.tsv",sep="\t")

#select top 3 genes out of 5 or 10 which are common in between both test as a candidate gene and compare with other exons from same gene.
plotSplice(ex, coef=ncol(ex), geneid="MSTRG.2185", genecol="GeneID", FDR = 0.05)
plotSplice(ex, coef=ncol(ex), geneid="MSTRG.10818", genecol="GeneID", FDR = 0.05)
ballgown::geneIDs(bg)[16745]

## Plotting setup
tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow', 'cyan')
palette(tropical)

## Plotting gene abundance distribution after filtering low abundent genes with a variance across the samples of less than one
fpkm <- texpr(bg_filt, meas='FPKM')
fpkm <- log2(fpkm+1)
colnames(fpkm) <- c("NoSalt_r1", "NoSalt_r2","NoSalt_r3","Salt_r1", "Salt_r2", "Salt_r3", "Bac_r1" , "Bac_r2" , "Bac_r3" , "Metaba_r1" , "Metaba_r2" , "Metaba_r3")
##boxplot(fpkm, col=as.numeric(pheno_data$type), las=3,ylab='log2(FPKM+1)')






# Display the gene name for a single row of data 
ballgown::geneIDs(bg)[9622]
ballgown::transcriptIDs(bg)[9622]
fpkm <- texpr(bg_filt, meas='FPKM')
fpkm <- log2(fpkm+1)
# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
plot(fpkm[9622,] ~ pheno_data$type, border=c(2,3), main=paste(ballgown::geneIDs(bg)[9622],' : ', ballgown::transcriptIDs(bg)[9622]),pch=19, xlab="Type", ylab="log2(fpkm+1)")
points(fpkm[9622,] ~ jitter(as.numeric(pheno_data$type)), col=as.numeric(pheno_data$type))

# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[9622], bg, main=c('Gene in sample Salt_rep1'), sample=c('Salt_rep1'))
plotTranscripts(ballgown::geneIDs(bg)[9622], bg, main=c('Gene in sample NoSalt_rep1'), sample=c('NoSalt_rep1'))

plotTranscripts(ballgown::geneIDs(bg)[5980], gown=bg, samples='Salt_rep1', 
                meas='FPKM', colorby='transcript', 
                main='transcripts from gene 9622: Salt_rep1, FPKM')
plotTranscripts(ballgown::geneIDs(bg)[5980], bg, 
                samples=c('Salt_rep1', 'Salt_rep2','Salt_rep3','NoSalt_rep1', 'NoSalt_rep2', 'NoSalt_rep3'), 
                meas='FPKM', colorby='transcript')

plotMeans(ballgown::geneIDs(bg)[5980],bg,groupvar="type", meas='FPKM', colorby='transcript')
plotMeans(ballgown::geneIDs(bg)[16745],bg,groupvar="type", meas='FPKM', colorby='transcript')


plotMeans(gene='MSTRG.2185',bg,groupvar="type", meas='FPKM', colorby='transcript')
clusterTranscripts(gene='MSTRG.2185', gown=bg, k=2, method='kmeans')
plotLatentTranscripts(gene='MSTRG.2185', gown=bg, k=2, method='kmeans', returncluster=FALSE)
agg = collapseTranscripts(gene='MSTRG.2185', gown=bg, k=2, method='kmeans')
stattest(gowntable=agg$tab, pData=pData(bg), feature='transcript_cluster', 
         covariate='type', libadjust=FALSE)
#for multiple gene list visualization
lapply(list_of_genes, plotTranscripts, gown=bg, samples = c('OLD_PLUS', 'OLD_MINUS',
                                                            'YOUNG_PLUS', 'YOUNG_MINUS') )