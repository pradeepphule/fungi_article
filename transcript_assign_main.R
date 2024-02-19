library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(supraHex)
library(limma)
#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))
#Set working directory where results files exist
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown" 
setwd(working_dir)
pheno_data = read.csv("all_path.csv")
descriptions = read.delim("transcripts_with_descriptions.txt") # not needed
trans_annotation = read.delim("blast2go_annotation_descriptions_20170509_1442.txt") 
trans_blastgotable = read.delim("blast2go_go_table_20170509_1402.txt")#complete annotation from blast2go
trans_kegg = read.delim("blast2go_kegg_20170509_1458.txt")# Kegg specific annotation
trans_annotation$Enzymes = trans_blastgotable$Enzyme.Names[match(trans_annotation$SeqName,trans_blastgotable$SeqName)]
trans_annotation$Pathway = trans_kegg$Pathway[match(trans_annotation$SeqName,trans_kegg$Seq)]
head(trans_annotation)
tmap_stringtie = read.delim("mapped_file.tmap")
head(tmap_stringtie)
trans_annotation <- merge(trans_annotation, tmap_stringtie,by.x=c("SeqName"),by.y=c("transcript_id"),all.x=TRUE)
head(trans_annotation)
########################################################################################################################
bg_all = ballgown(samples=as.vector(pheno_data$path),pData=pheno_data, meas='all')
bg_all # ballgown instance with 17122 transcripts and 12 samples
bg_all_table = texpr(bg_all, 'all')
head(bg_all_table)
bg_all_filt = subset(bg_all,"rowVars(texpr(bg_all)) > 1", genomesubset=TRUE)#Low abundance genes were filtered by removing transcripts with a variance across the sample less than one. 
bg_all_filt # ballgown instance with 11520 transcripts and 12 samples

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg_all_filt, 'all')
head(bg_table)
bg_gene_names = unique(bg_table[, 9:10])
head(bg_gene_names)
dim(bg_gene_names)#7842 genes with    2 columnns
# start from here after running ballgown_main scriptPull the gene_expression data frame from the ballgown object
transcript_expression = as.data.frame(texpr(bg_all_filt))
head(transcript_expression)
colnames(transcript_expression) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
head(transcript_expression)
dim(transcript_expression)
min_nonzero=1
data_columns=c(1:12)
data_colors=c("hotpink1", "hotpink2", "hotpink3","tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3","cyan1","cyan2","cyan3")

boxplot(log2(transcript_expression[,data_columns]+min_nonzero), col=data_colors, las=3, ylab="log2(FPKM+1)", main="Transcript level abundance distribution ")
# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg_all_filt)$t2g
transcript_gene_table$transcriptname = bg_table$t_name[match(transcript_gene_table$t_id,bg_table$t_id)]
head(transcript_gene_table)
dim(transcript_gene_table)
transcript_gene_table = merge(transcript_gene_table,trans_annotation,by.x=c("transcriptname"),by.y=c("SeqName"))
head(transcript_gene_table)
dim(transcript_gene_table)
unique_gene_id = unique(transcript_gene_table$g_id)
head(unique_gene_id)
ordered_gene_id <- with(transcript_gene_table, transcript_gene_table[order(g_id, transcriptname, -t_id),])
head(ordered_gene_id)
dim(ordered_gene_id)
gene_annotation  = ordered_gene_id[!duplicated(ordered_gene_id$g_id), ]
head(gene_annotation)
dim(gene_annotation)
dim(trans_annotation)
dim(transcript_gene_table)
dim(bg_table)
write.table(transcript_gene_table,"transcript_annotation_filtered.tsv",sep="\t")#11520 transcript annotation
write.table(gene_annotation,"gene_annotation_all.tsv",sep="\t")#7842 gene annotation
write.table(bg_table,"bg_table_all_filtered.tsv",sep="\t")# fpkm expression of 11520 transcripts.
write.table(trans_annotation,"transcript_annotation_all.tsv",sep="\t")#17122 transcript with 5 columns
write.table(bg_all_table,"bg_table_all.tsv",sep="\t")# fpkm expression of 17122 transcripts.
transcript_complete_annotation = merge(trans_annotation,bg_all_table,by.x=c("SeqName"),by.y=c("t_name"))
head(transcript_complete_annotation)
write.table(transcript_complete_annotation,"transcript_complete_annotation.tsv",sep="\t")# fpkm expression of 17122 transcripts with annotation
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
###################################################################################################################################################################
###################################################################################################################################################################

bg1 = subset(bg_all_filt, "pheno_data$type =='Salt' | pheno_data$type =='NoSalt'",genomesubset = FALSE)
#bg1 = ballgown(samples=as.vector(pheno_data1$path),pData=pheno_data1, meas='all')
bg1
#bg1_filt = subset (bg1,"rowVars(texpr(bg1)) > 1", genomesubset=TRUE)
#bg1_filt
bg1_table = texpr(bg1, 'all')
head(bg1_table)
dim(bg1_table)
bg1_gene_names = unique(bg1_table[, 9:10])
rt1 = stattest(bg1, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
rt1 <- data.frame(rt1,log_fc = log2(rt1$fc),Log10_Pvalue = log10(rt1$pval))
#rt1 = merge(rt1,bg1_table,by.x=c("id"),by.y=c("t_id"))
rt1$feature = bg1_table$t_name[match(rt1$id,bg1_table$t_id)]
rt1$id <- NULL
rt1 = merge(rt1,trans_annotation,by.x=c("feature"),by.y=c("SeqName"))
dim(rt1)
head(rt1)
sig1 = subset(rt1,rt1$pval<0.05)
sig1 = arrange(sig1,pval)
head(sig1)
dim(sig1)
## get the p-values and calculate the scores thereupon
sig1_pval <- sig1$pval
## look at the distribution of p-values
hist(sig1_pval, breaks=50, xlab="Pval", main="Histogram of significant pvalues SaltVsNosalt", col="skyblue")
write.table(sig1,"Salt_vs_NoSalt_transcript_results_sig_arrangedbypval.tsv",sep="\t")
# Plot #4 - View the distribution of differential expression values as a histogram
#Display only those that are significant according to Ballgown fold change greater than 4
hist(sig1[,"log_fc"], breaks=50, col="seagreen", xlab="log2(Fold change) SaltVsNoSalt", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 2", lwd=2, lty=2)

sig1_logfc2 = subset(sig1,abs(sig1[,"log_fc"] )>= 2)
sig1_logfc2 = arrange(sig1_logfc2,desc(log_fc))
write.table(sig1_logfc2,"Salt_vs_NoSalt_transcript_results_sig_arrangedbylogfc2.tsv",sep="\t")#use it for venn diagram.
# to get top level result just order by p value no need of order by fc and then q value , it can be done in excel
sig1 = arrange(sig1_logfc2,pval)
head(sig1_logfc2)
dim(sig1_logfc2)
#select top 25 entries for comparision at later stage. this output can be used for top level comparision
write.table(sig1_logfc2,"Salt_vs_NoSalt_transcript_results_sig_arrangedbylogfc2_pval.tsv",sep="\t")



bg2 = subset(bg_all_filt, "pheno_data$type =='Salt' | pheno_data$type =='Bac'",genomesubset = FALSE)
#bg2 = ballgown(samples=as.vector(pheno_data2$path),pData=pheno_data2, meas='all')
bg2
#bg2_filt = subset (bg2,"rowVars(texpr(bg2)) > 1", genomesubset=TRUE)
#bg2_filt
bg2_table = texpr(bg2, 'all')
bg2_gene_names = unique(bg2_table[, 9:10])
rt2 = stattest(bg2, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
rt2 <- data.frame(rt2,log_fc = log2(rt2$fc),Log10_Pvalue = log10(rt2$pval))
#rt2 = merge(rt2,bg2_table,by.x=c("id"),by.y=c("t_id"))
rt2$feature = bg2_table$t_name[match(rt2$id,bg2_table$t_id)]
rt2$id <- NULL
rt2 = merge(rt2,trans_annotation,by.x=c("feature"),by.y=c("SeqName"))
dim(rt2)
head(rt2)
sig2 = subset(rt2,rt2$pval<0.05)
sig2 = arrange(sig2,pval)
head(sig2)
dim(sig2)
## get the p-values and calculate the scores thereupon
sig2_pval <- sig2$pval
## look at the distribution of p-values
hist(sig2_pval, breaks=50, xlab="Pval", main="Histogram of significant pvalues SaltVsBac", col="skyblue")
write.table(sig2,"Salt_vs_Bac_transcript_results_sig_arrangedbypval.tsv",sep="\t")

hist(sig2[,"log_fc"], breaks=50, col="seagreen", xlab="log2(Fold change) SaltVsBac", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 2", lwd=2, lty=2)

sig2_logfc2= subset(sig2,abs(sig2[,"log_fc"] )>= 2)
sig2_logfc2 = arrange(sig2_logfc2,desc(log_fc))
write.table(sig2_logfc2,"Salt_vs_Bac_transcript_results_sig_arrangedbylogfc2.tsv",sep="\t")
sig2 = arrange(sig2_logfc2,pval)
head(sig2_logfc2)
dim(sig2_logfc2)
write.table(sig2_logfc2,"Salt_vs_Bac_transcript_results_sig_arrangedbylogfc2_pval.tsv",sep="\t")


bg3 = subset(bg_all_filt, "pheno_data$type =='Salt' | pheno_data$type =='Metabact'",genomesubset = FALSE)
#bg3 = ballgown(samples=as.vector(pheno_data3$path),pData=pheno_data3, meas='all')
bg3
#bg3_filt = subset (bg3,"rowVars(texpr(bg3)) > 1", genomesubset=TRUE)
#bg3_filt
bg3_table = texpr(bg3, 'all')
bg3_gene_names = unique(bg3_table[, 9:10])
rt3 = stattest(bg3, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
rt3 <- data.frame(rt3,log_fc = log2(rt3$fc),Log10_Pvalue = log10(rt3$pval))
#rt3 = merge(rt3,bg3_table,by.x=c("id"),by.y=c("t_id"))
rt3$feature = bg3_table$t_name[match(rt3$id,bg3_table$t_id)]
rt3$id <- NULL
rt3 = merge(rt3,trans_annotation,by.x=c("feature"),by.y=c("SeqName"))
dim(rt3)
head(rt3)
sig3 = subset(rt3,rt3$pval<0.05)
sig3 = arrange(sig3,pval)
head(sig3)
dim(sig3)
## get the p-values and calculate the scores thereupon
sig3_pval <- sig3$pval
## look at the distribution of p-values
hist(sig3_pval, breaks=50, xlab="Pval", main="Histogram of significant pvalues SaltVsMetabact", col="skyblue")
write.table(sig3,"Salt_vs_Metabact_transcript_results_sig_arrangedbypval.tsv",sep="\t")

hist(sig3[,"log_fc"], breaks=50, col="seagreen", xlab="log2(Fold change) SaltVsMetabact", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 2", lwd=2, lty=2)

sig3_logfc2 = subset(sig3,abs(sig3[,"log_fc"] )>= 2)
sig3_logfc2 = arrange(sig3_logfc2,desc(log_fc))
write.table(sig3_logfc2,"Salt_vs_Metabact_transcript_results_sig_arrangedbylogfc2.tsv",sep="\t")
sig3 = arrange(sig3_logfc2,pval)
head(sig3_logfc2)
dim(sig3_logfc2)
write.table(sig3_logfc2,"Salt_vs_Metabact_transcript_results_sig_arrangedbylogfc2_pval.tsv",sep="\t")

bg4 = subset(bg_all_filt, "pheno_data$type =='Bac' | pheno_data$type =='Metabact'",genomesubset = FALSE)
#bg4 = ballgown(samples=as.vector(pheno_data4$path),pData=pheno_data4, meas='all')
bg4
#bg4_filt = subset (bg4,"rowVars(texpr(bg4)) > 1", genomesubset=TRUE)
#bg4_filt
bg4_table = texpr(bg4, 'all')
bg4_gene_names = unique(bg4_table[, 9:10])
rt4 = stattest(bg4, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
rt4 <- data.frame(rt4,log_fc = log2(rt4$fc),Log10_Pvalue = log10(rt4$pval))
#rt4 = merge(rt4,bg4_table,by.x=c("id"),by.y=c("t_id"))
rt4$feature = bg4_table$t_name[match(rt4$id,bg4_table$t_id)]
rt4$id <- NULL
rt4 = merge(rt4,trans_annotation,by.x=c("feature"),by.y=c("SeqName"))
dim(rt4)
head(rt4)
sig4 = subset(rt4,rt4$pval<0.05)
sig4 = arrange(sig4,pval)
head(sig4)
dim(sig4)
## get the p-values and calculate the scores thereupon
sig4_pval <- sig4$pval
## look at the distribution of p-values
hist(sig4_pval, breaks=50, xlab="Pval", main="Histogram of significant pvalues BacVsMetabact", col="skyblue")
write.table(sig4,"Bac_vs_Metabact_transcript_results_sig_arrangedbypval.tsv",sep="\t")

hist(sig4[,"log_fc"], breaks=50, col="seagreen", xlab="log2(Fold change) BacVsMetabact", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 2", lwd=2, lty=2)

sig4_logfc2 = subset(sig4,abs(sig4[,"log_fc"] )>= 2)
sig4_logfc2 = arrange(sig4_logfc2,desc(log_fc))
write.table(sig4_logfc2,"Bac_vs_Metabact_transcript_results_sig_arrangedbylogfc2.tsv",sep="\t")
sig4_logfc2 = arrange(sig4_logfc2,pval)
head(sig4_logfc2)
dim(sig4_logfc2)
write.table(sig4_logfc2,"Bac_vs_Metabact_transcript_results_sig_arrangedbylogfc2_pval.tsv",sep="\t")


#### Plot #5 - Display the grand expression values from carcinoma and normal and mark those that are significantly differentially expressed
transcript_expression[,"Salt"]=apply(transcript_expression[,c(1:3)], 1, mean)
transcript_expression[,"NoSalt"]=apply(transcript_expression[,c(4:6)], 1, mean)
transcript_expression[,"Bac"]=apply(transcript_expression[,c(7:9)], 1, mean)
transcript_expression[,"Metabact"]=apply(transcript_expression[,c(10:12)], 1, mean)
head(transcript_expression)
a=log2(transcript_expression[,"Salt"]+min_nonzero)
b=log2(transcript_expression[,"NoSalt"]+min_nonzero)
c=log2(transcript_expression[,"Bac"]+min_nonzero)
d=log2(transcript_expression[,"Metabact"]+min_nonzero)

sig1_index=which(rt1$pval<0.05)
plot(x=a, y=b, pch=16, cex=0.25, xlab="Salt FPKM (log2)", ylab="NoSalt FPKM (log2)", main="SaltVsNOSalt FPKMs")
abline(a=0, b=1)
asig=a[sig1_index]
bsig=b[sig1_index]
points(x=asig, y=bsig, col="magenta", pch=16, cex=0.25)
legend("topleft", "Significant", col="magenta", pch=16)

sig2_index=which(rt2$pval<0.05)
plot(x=a, y=c, pch=16, cex=0.25, xlab="Salt FPKM (log2)", ylab="Bac FPKM (log2)", main="SaltVsBac FPKMs")
abline(a=0, b=1)
asig=a[sig2_index]
csig=c[sig2_index]
points(x=asig, y=csig, col="magenta", pch=16, cex=0.25)
legend("topleft", "Significant", col="magenta", pch=16)

sig3_index=which(rt3$pval<0.05)
plot(x=a, y=d, pch=16, cex=0.25, xlab="Salt FPKM (log2)", ylab="Metabact FPKM (log2)", main="SaltVsMetabact FPKMs")
abline(a=0, b=1)
asig=a[sig3_index]
dsig=d[sig3_index]
points(x=asig, y=dsig, col="magenta", pch=16, cex=0.25)
legend("topleft", "Significant", col="magenta", pch=16)

sig4_index=which(rt4$pval<0.05)
plot(x=c, y=d, pch=16, cex=0.25, xlab="Bac FPKM (log2)", ylab="Metabact FPKM (log2)", main="BacVsMetabact FPKMs")
abline(a=0, b=1)
csig=c[sig1_index]
dsig=d[sig1_index]
points(x=csig, y=dsig, col="magenta", pch=16, cex=0.25)
legend("topleft", "Significant", col="magenta", pch=16)


#############################################################################################################
#################################################################################################################
#################################################################################################################

# combine all diff expressed transcripts into one table
df <- rbind(sig1,sig2,sig3,sig4)
df_logfc2 <- rbind(sig1_logfc2,sig2_logfc2,sig3_logfc2,sig4_logfc2)
head(df_logfc2)
nrow(df_logfc2)
# use this command to select unique transcript id from the diff. expressed transcripts.
u <- unique(df[,1])
u_df_logfc2 <- unique(df_logfc2[,1])

# retrive transcript level fpkm values from transcript ids. and make matrix of these ids.

indices <- match( unique(df[,1]), texpr(bg_all, 'all')$t_name)
t_n <- texpr(bg_all, 'all')$t_name[indices]

indices_logfc2 <- match( unique(df_logfc2[,1]), texpr(bg_all, 'all')$t_name)
t_n_logfc2 <- texpr(bg_all, 'all')$t_name[indices_logfc2]

#make a table of all transcriptid and rpkm values from bg_all table by matching id with table
  #it should be 1 column of t_name and all 12 sample fpkm values
  #convert t_name into rownames 
  #convert colnames into rownames
#rownames(df) <- df$ID
#i = row.names(gene_expression) == "transcipt_id_object"
#gene_expression[i,]
#name_df = read.table('gene_file_name')
#name_column = name_df[,1]
bg_sig = subset(bg_all, "t_name %in% t_n")
#bg_sig = subset(bg_all, "t_name %in% t_n_logfc2")
sig_trans_all <- as.data.frame(texpr(bg_sig, meas='FPKM'))
rownames(sig_trans_all) <- t_n
#rownames(sig_trans_all) <- t_n_logfc2
head(sig_trans_all)
dim(sig_trans_all)
#check above transcript manually in table
bg_sig_table = texpr(bg_sig, 'all')
head(bg_sig_table)
colnames(sig_trans_all) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
data_columns=c(1:12)
min_nonzero=1
#to normalize with log2-transform and then mean-center the data for creating heatmap based on RPKM.
sig_trans_all = log2(sig_trans_all[,data_columns]+min_nonzero)
head(sig_trans_all)
# transform data by row/gene centering
data_sigtrans <- sig_trans_all - matrix(rep(apply(sig_trans_all,1,mean),ncol(sig_trans_all)),ncol=ncol(sig_trans_all))
head(data_sigtrans)
# (II) Train the supra-hexagonal map with input data only
sMap <- sPipeline(data_sigtrans)
visHexMulComp(sMap, title.rotate=10)
sWriteData(sMap, data_sigtrans, filename="Output_sMap_Bg_data_sigtransall.txt")
# (IV) Perform partitioning operation on the map to obtain continuous clusters (i.e. gene meta-clusters) as they are different from gene clusters in an individual map node
sBase <- sDmatCluster(sMap)
visDmatCluster(sMap, sBase)
sWriteData(sMap, data_sigtrans, sBase, filename="Output_sBase_Bg_data_sigtransall.txt")
# prepare colors for the column sidebar
# color for stages (S9-S14)
stages <- sub("_.*","",colnames( data_sigtrans))
lvs <- unique(stages)
lvs_color <- visColormap(colormap="jet")(length(lvs))
col_stages <- sapply(stages, function(x) lvs_color[x==lvs])
# color for replicates (R1-R3)
replicates <- sub(".*_","",colnames(data_sigtrans))
lvs <- unique(replicates)
lvs_color <- visColormap(colormap="gray-black")(length(lvs))
col_replicates <- sapply(replicates, function(x) lvs_color[x==lvs])
# combine both color vectors
ColSideColors <- cbind(col_stages,col_replicates)
colnames(ColSideColors) <- c("conditions","Replicates")

#pdf(file="Bg_data_sigtransall.pdf")
output <- visDmatHeatmap(sMap,  data_sigtrans, sBase, base.legend.location="bottomleft", reorderRow="hclust", ColSideColors=ColSideColors, KeyValueName="log2(Ratio)", ColSideLabelLocation="right", labRow=NA)
output
dev.off()
# (V) Reorder the sample-specific components of the map to delineate relationships between samples
sReorder <- sCompReorder(data_sigtrans, metric="euclidean")
# svg("AAPL.svg",width=14,height=7)
visCompReorder(sMap, sReorder, title.rotate=15)
dev.off()

################################################################################
###################plots and splicing analysis##################################
################################################################################
library(limma)
limmaUsersGuide()
?topSplice()
#fit a linear model from exon level read count data expr(bg1)gives read count data which is fitted inside linear model
#here exon level splicing analysis is done by read count data.
#You could use limma instead: fit linear models on each row of the exon count table (eexpr(bg) in ballgown), then use the diffSplice function to test for differential splicing.
fit1 <- lmFit(eexpr(bg1), c(0,0,0,1,1,1))
dim(fit1)
#exon/transcript mapping---all (unique/first by strart position) exon id are matched with duplicated/multiple exons of a transcript and transcript ids are selected.
transcript_id_by_exon1 = indexes(bg1)$e2t$t_id[match(unique(indexes(bg1)$e2t$e_id), indexes(bg1)$e2t$e_id)]
#transcript/gene mapping---all transcript ids are seleted.including duplicated/multiple for sigle gene.
transcript_id_by_gene1 = indexes(bg1)$t2g$t_id #transcript/gene mapping
# gene id is retrived by matching unique exon id containing transcript with transcript id of the gene. 
geneID1 = indexes(bg1)$t2g$g_id[match(transcript_id_by_exon1, transcript_id_by_gene1)]#gene id s are calc. by matching t_id of exon to t_id of gene
# collect results from diffsplice function was used on count based linear model with geneID1
ex1 <- diffSplice(fit1, geneID1)
#top 20 geneid, no.of exon , pvalue, fdr splicing results
#write.table(ex1,"Salt_vs_NoSalt_diffSplice_counts.tsv",sep="\t")

##for getting top level splicing from result ////use simes test which calculate t test for all exons of the gene and returns its p value.
#The simes tests for any differences in exon usage between experimental conditions.
#"t" gives t-tests for each exon. "simes" gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
##if exon no. are more from same gene then there is more possibility of differentially spliceing from same gene in differnet conditions.
#topsplice_simes1<-topSplice(ex1, test="simes", n=20,FDR=1)
#topsplice_simes1 = merge(topsplice_simes1,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
#topsplice_simes1
#top 20 geneid, no.of exon , pvalue, fdr splicing results
#write.table(topsplice_simes1,"Salt_vs_NoSalt_topsplice_simes_counts.tsv",sep="\t")
#"F" gives F-tests for each gene.F	moderated F-statistic (if level="gene")
topsplice_ftest1<-topSplice(ex1, test="F", n=100,FDR=1)
topsplice_ftest1 = merge(topsplice_ftest1,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
#topsplice_ftest1$t_id <- NULL
#topsplice_ftest1$transcriptname <- NULL
topsplice_ftest1
#top 20 geneid, no.of exon ,f test,  pvalue, fdr splicing results
write.table(topsplice_ftest1,"Salt_vs_NoSalt_topsplice_ftest_counts.tsv",sep="\t")
topsplice_ftest1_200<-topSplice(ex1, test="F", n=200,FDR=1)
topsplice_ftest1_200 = merge(topsplice_ftest1_200,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
write.table(topsplice_ftest1_200 ,"Salt_vs_NoSalt_topsplice_ftest_counts_200.tsv",sep="\t")

#fit a linear model from exon level read count data
fit2 <- lmFit(eexpr(bg2), c(0,0,0,1,1,1))
dim(fit2)
transcript_id_by_exon2 = indexes(bg2)$e2t$t_id[match(unique(indexes(bg2)$e2t$e_id), indexes(bg2)$e2t$e_id)]
transcript_id_by_gene2 = indexes(bg2)$t2g$t_id #transcript/gene mapping
geneID2 = indexes(bg2)$t2g$g_id[match(transcript_id_by_exon2, transcript_id_by_gene2)]
ex2 <- diffSplice(fit2, geneID2)
#topsplice_simes2<-topSplice(ex2, test="simes", n=20,FDR=1)
#topsplice_simes2 = merge(topsplice_simes2,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
#topsplice_simes2
#write.table(ex2,"SaltVsBac_diffSplice_counts.tsv",sep="\t")
topsplice_ftest2<-topSplice(ex2, test="F", n=100,FDR=1)
topsplice_ftest2 = merge(topsplice_ftest2,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
topsplice_ftest2
write.table(topsplice_ftest2,"SaltVsBac_topsplice_ftest_counts.tsv",sep="\t")
topsplice_ftest2_200<-topSplice(ex2, test="F", n=200,FDR=1)
topsplice_ftest2_200 = merge(topsplice_ftest2_200,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
write.table(topsplice_ftest2_200 ,"SaltVsBac_topsplice_ftest_counts_200.tsv",sep="\t")

fit3 <- lmFit(eexpr(bg3), c(0,0,0,1,1,1))
dim(fit3)
transcript_id_by_exon3 = indexes(bg3)$e2t$t_id[match(unique(indexes(bg3)$e2t$e_id), indexes(bg3)$e2t$e_id)]
transcript_id_by_gene3 = indexes(bg3)$t2g$t_id #transcript/gene mapping
geneID3 = indexes(bg3)$t2g$g_id[match(transcript_id_by_exon3, transcript_id_by_gene3)]#fit a linear model from exon level read count data
ex3 <- diffSplice(fit3, geneID3)
#topsplice_simes3<-topSplice(ex3, test="simes", n=20,FDR=1)
#topsplice_simes3 = merge(topsplice_simes3,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
#topsplice_simes3
#write.table(ex3,"SaltVsMetabact_diffSplice_counts.tsv",sep="\t")
topsplice_ftest3<-topSplice(ex3, test="F", n=100,FDR=1)
topsplice_ftest3 = merge(topsplice_ftest3,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
topsplice_ftest3
write.table(topsplice_ftest3,"SaltVsMetabact_topsplice_ftest_counts.tsv",sep="\t")
topsplice_ftest3_200<-topSplice(ex3, test="F", n=200,FDR=1)
topsplice_ftest3_200 = merge(topsplice_ftest3_200,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
write.table(topsplice_ftest3_200 ,"SaltVsMetabact_topsplice_ftest_counts_200.tsv",sep="\t")

fit4 <- lmFit(eexpr(bg4), c(0,0,0,1,1,1))
dim(fit4)
transcript_id_by_exon4 = indexes(bg4)$e2t$t_id[match(unique(indexes(bg4)$e2t$e_id), indexes(bg4)$e2t$e_id)]
transcript_id_by_gene4 = indexes(bg4)$t2g$t_id #transcript/gene mapping
geneID4 = indexes(bg4)$t2g$g_id[match(transcript_id_by_exon4, transcript_id_by_gene4)]
ex4 <- diffSplice(fit4, geneID4)
#topsplice_simes4<-topSplice(ex4, test="simes", n=20,FDR=1)
#topsplice_simes4 = merge(topsplice_simes4,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
#topsplice_simes4
#write.table(ex4,"BacVsMetabact_diffSplice_counts.tsv",sep="\t")
topsplice_ftest4<-topSplice(ex4, test="F", n=100,FDR=1)
topsplice_ftest4 = merge(topsplice_ftest4,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
topsplice_ftest4
write.table(topsplice_ftest4,"BacVsMetabact_topsplice_ftest_counts.tsv",sep="\t")
topsplice_ftest4_200<-topSplice(ex4, test="F", n=200,FDR=1)
topsplice_ftest4_200 = merge(topsplice_ftest4_200,gene_annotation,by.x=c("GeneID"),by.y=c("g_id"))
write.table(topsplice_ftest4_200 ,"BacVsMetabact_topsplice_ftest_counts_200.tsv",sep="\t")


working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/results/result_pval/results_logfc2/results_logfc2_arrangedby_pval/splicing/SaltVsNoSalt_topsplice_selected_graphs" 
setwd(working_dir)
topsplice_ftest1_200_selected<-read.table("Salt_vs_NoSalt_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
topsplice_ftest1_200_selected_GeneID<- as.matrix(topsplice_ftest1_200_selected$GeneID)
lapply(topsplice_ftest1_200_selected_GeneID, plotMeans,gown=bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)

working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/results/result_pval/results_logfc2/results_logfc2_arrangedby_pval/splicing/SaltVsBac_topsplice_selected_graphs" 
setwd(working_dir)
topsplice_ftest2_200_selected<-read.table("SaltVsBac_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
topsplice_ftest2_200_selected_GeneID<- as.matrix(topsplice_ftest2_200_selected$GeneID)
lapply(topsplice_ftest2_200_selected_GeneID, plotMeans,gown=bg2,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)


working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/results/result_pval/results_logfc2/results_logfc2_arrangedby_pval/splicing/SaltVsMetabact_topsplice_selected_graphs" 
setwd(working_dir)
#plotSplice(ex4, geneid="Fpp1-1_1088", genecol="GeneID")
topsplice_ftest3_200_selected<-read.table("SaltVsMetabact_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
topsplice_ftest3_200_selected_GeneID<- as.matrix(topsplice_ftest3_200_selected$GeneID)
lapply(topsplice_ftest3_200_selected_GeneID, plotMeans,gown=bg3,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)




working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/ballgown/results/result_pval/results_logfc2/results_logfc2_arrangedby_pval/splicing/BacVsMetabact_topsplice_selected_graphs" 
setwd(working_dir)
topsplice_ftest4_200_selected<-read.table("BacVsMetabact_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
topsplice_ftest4_200_selected_GeneID<- as.matrix(topsplice_ftest4_200_selected$GeneID)
lapply(topsplice_ftest4_200_selected_GeneID, plotMeans,gown=bg4,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)




#topsplice_ftest2_200_selected<-read.table("SaltVsBac_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
#topsplice_ftest3_200_selected<-read.table("SaltVsMetabact_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")
#topsplice_ftest4_200_selected<-read.table("BacVsMetabact_topsplice_ftest_counts_200_selected.tsv",header = TRUE, sep="\t")


#topsplice_ftest2_200_selected_GeneID<- as.matrix(topsplice_ftest2_200_selected$GeneID)
#topsplice_ftest3_200_selected_GeneID<- as.matrix(topsplice_ftest3_200_selected$GeneID)
#topsplice_ftest4_200_selected_GeneID<- as.matrix(topsplice_ftest4_200_selected$GeneID)


#lapply(topsplice_ftest2_200_selected_GeneID, plotMeans,gown=bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)
#lapply(topsplice_ftest3_200_selected_GeneID, plotMeans,gown=bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)
#lapply(topsplice_ftest4_200_selected_GeneID, plotMeans,gown=bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)




##########lapply(topsplice_ftest4$GeneID, plotTranscripts, gown=bg4,samples = c('Bac_rep1', 'Bac_rep1',
                                                                      'Bac_rep1'))
###########lapply(topsplice_ftest4$GeneID, plotMeans, gown=bg4,groupvar="type", meas='FPKM', colorby='transcript')
#lapply(list_of_genes, plotTranscripts, gown=bg, samples = c('OLD_PLUS', 'OLD_MINUS','YOUNG_PLUS', 'YOUNG_MINUS') )

#plotSplice(ex1,coef=ncol(ex1), geneid="Fpp1-1_6509", genecol="GeneID",FDR = 0.05)
########plotMeans(gene='Fpp1-1_1133',bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)
########plotMeans(gene='Fpp1-1_6509',bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)
########plotMeans(gene='Fpp1-1_6510',bg1,groupvar="type",meas='rcount', colorby='exon',labelTranscripts = TRUE)
#plotMeans(gene='Fpp1-1_736',bg4,groupvar="type", meas='FPKM', colorby='transcript')
#plotMeans(gene='Fpp1-1_1088',bg_all_filt,groupvar="type", meas='FPKM', colorby='transcript')
#plotMeans(gene='Fpp1-1_736',bg_all_filt,groupvar="type", meas='FPKM', colorby='transcript')
#plotMeans(gene='Fpp1-1_2185',bg_all_filt,groupvar="type", meas='FPKM', colorby='transcript')
# transcript clustering
#clusterTranscripts(gene='Fpp1-1_2185', gown=bg_all_filt, k=3, method='kmeans')
#plotLatentTranscripts(gene='Fpp1-1_2185', gown=bg_all_filt, k=3, method='kmeans', returncluster=FALSE)
#to calculate p value for significance of transcript clusters from same gene , --needed when you want to do validation.
#agg = collapseTranscripts(gene='Fpp1-1_2185', gown=bg_all_filt, k=3, method='kmeans')
#stattest(gowntable=agg$tab, pData=pData(bg_all_filt), feature='transcript_cluster', 
#        covariate='type', libadjust=FALSE)

#####similarly all secreated, effector and cystein rich , transcription factors, cazy, transporter, cytochrome genes protein encoding genes can be found.
####find out only interproscan id for transcripts. match with transcriptdb and finalize transcrption factor proteins.
