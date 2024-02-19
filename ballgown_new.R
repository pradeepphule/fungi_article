library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(supraHex)
#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))

#Set working directory where results files exist
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown" 
setwd(working_dir)
#pdf(file="output.pdf")

pheno_data = read.csv("all_path.csv")

# Load ballgown data structure and save it to a variable "bg"
bg_all = ballgown(samples=as.vector(pheno_data$path),pData=pheno_data, meas='all')


# Display a description of this object
bg_all
# Save the ballgown object to a file for later use
save(bg_all, file='bg_all.rda')

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_all_filt = subset (bg_all,"rowVars(texpr(bg_all)) > 1", genomesubset=TRUE)

bg_all_filt
# Load gene names for lookup later in the tutorial
bg_table = texpr(bg_all_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])


# start from here after running ballgown_main scriptPull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg_all_filt))
colnames(gene_expression) <- c("Salt.r1","Salt.r2","Salt.r3","NOSalt.r1","NOSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
head(gene_expression)
dim(gene_expression)
min_nonzero=1
data_columns=c(1:12)
data_colors=c("hotpink1", "hotpink2", "hotpink3","tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3","cyan1","cyan2","cyan3")

boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, las=3, ylab="log2(FPKM+1)", main="Gene level abundance distribution ")
#to normalize with log2-transform and then mean-center the data for creating heatmap based on RPKM.
gene_expression = log2(gene_expression[,data_columns]+min_nonzero)
head(gene_expression)
# transform data by row/gene centering
data_gene <- gene_expression - matrix(rep(apply(gene_expression,1,mean),ncol(gene_expression)),ncol=ncol(gene_expression))
head(data_gene)
# (II) Train the supra-hexagonal map with input data only
sMap <- sPipeline(data_gene)
visHexMulComp(sMap, title.rotate=15)
dev.off()
sWriteData(sMap, data_gene, filename="Output_sMap_Bg_all.txt")
# visHexMapping(sMap, mappingType="indexes")
# visHexMapping(sMap, mappingType="hits")
# visHexMapping(sMap, mappingType="dist")
# visHexPattern(sMap, plotType="lines")
# visHexPattern(sMap, plotType="bars",colormap="rainbow",legend.cex=0.5)

# (IV) Perform partitioning operation on the map to obtain continuous clusters (i.e. gene meta-clusters) as they are different from gene clusters in an individual map node
sBase <- sDmatCluster(sMap)
visDmatCluster(sMap, sBase)
sWriteData(sMap, data_gene, sBase, filename="Output_sBase_Bg_all.txt")
# prepare colors for the column sidebar
# color for stages (S9-S14)
stages <- sub("_.*","",colnames(data_gene))
lvs <- unique(stages)
lvs_color <- visColormap(colormap="jet")(length(lvs))
col_stages <- sapply(stages, function(x) lvs_color[x==lvs])
# color for replicates (R1-R3)
replicates <- sub(".*_","",colnames(data_gene))
lvs <- unique(replicates)
lvs_color <- visColormap(colormap="gray-black")(length(lvs))
col_replicates <- sapply(replicates, function(x) lvs_color[x==lvs])
# combine both color vectors
ColSideColors <- cbind(col_stages,col_replicates)
colnames(ColSideColors) <- c("conditions","Replicates")
par(mar=rep(2,6))
pdf(file="heatmap_bg_all.pdf")
output <- visDmatHeatmap(sMap, data_gene, sBase, base.legend.location="bottomleft", reorderRow="hclust", ColSideColors=ColSideColors, KeyValueName="log2(Ratio)", ColSideLabelLocation="right", labRow=NA)
dev.off()
# (V) Reorder the sample-specific components of the map to delineate relationships between samples
sReorder <- sCompReorder(data_gene, metric="euclidean")
# svg("AAPL.svg",width=14,height=7)
visCompReorder(sMap, sReorder, title.rotate=15)
dev.off()
# exon level read count expression
# count_expression = as.data.frame(eexpr(bg_all_filt))
# head(count_expression)
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg_all_filt)$t2g
head(transcript_gene_table)

#### Plot #1 - the number of transcripts per gene.  
#Many genes will have only 1 transcript, some genes will have several transcripts
#Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol)
#Then use the 'hist' command to create a histogram of these counts
#How many genes have 1 transcript?  More than one transcript?  What is the maximum number of transcripts for a single gene?
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend_text
legend("topright", legend_text, lty=NULL)

#### Plot #2 - the distribution of transcript sizes as a histogram
#In this analysis we supplied StringTie with transcript models so the lengths will be those of known transcripts
#However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled
#If we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts

full_table <- texpr(bg_all_filt , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")