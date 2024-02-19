library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data from a file we saved in the current working directory
pheno_data = read.csv("/home/pradeep/RNA_HOME/fusarium/ballgown/rename/Salt_vs_NoSalt.csv")

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path),pData=pheno_data, meas='all')


# Display a description of this object
bg

# Load all attributes including gene name
#bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
bg_gene_info = unique(bg_table[, 1:10])
# Save the ballgown object to a file for later use
save(bg, file='bg.rda')

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_info,by.x=c("id"),by.y=c("gene_id"))
results_transcripts = merge(results_genes,bg_gene_info,by.x=c("id"),by.y=c("gene_id"))
# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_transcript_results.tsv",sep="\t")
write.table(results_genes,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_gene_results.tsv",sep="\t")

# get the mapping relationship between gene symbols and gene ids
# genename_ids <- dplyr::filter(results_transcripts, geneNames!=".") %>% dplyr::distinct(geneNames, geneIDs)
# join them by gene ids.
# results_genes1 <- dplyr::left_join(results_genes, genename_ids, by=c("id"="geneIDs"))



# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
bg_filt_allnames = unique(bg_filt_table)
# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
##results_genes = merge(results_genes,bg_filt_gene_names,by.x=c("id"),by.y=c("gene_id"))
## or same thing can be achieved by below code
## Add gene name
##results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt),
##                                 geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
results_transcripts = merge(results_transcripts,bg_filt_allnames,by.x=c("id"),by.y=c("gene_id"))
results_genes = merge(results_genes,bg_filt_allnames,by.x=c("id"),by.y=c("gene_id"))


## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)
# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts,"/home/pradeep/RNA_HOME/fusarium/ballgown/rename/Salt_vs_NoSalt_transcript_results_filtered.tsv",sep="\t")
write.table(results_genes,"/home/pradeep/RNA_HOME/fusarium/ballgown/rename/Salt_vs_NoSalt_gene_results_filtered.tsv",sep="\t")

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts,results_transcripts$pval<0.05)
sig_genes = subset(results_genes,results_genes$pval<0.05)
## Sort results from smallest p-value
sig_transcripts <- arrange(sig_transcripts, pval)
sig_genes <-  arrange(sig_genes, pval)

# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts,"/home/pradeep/RNA_HOME/fusarium/ballgown/rename/Salt_vs_NoSalt_transcript_results_sig.tsv",sep="\t")
write.table(sig_genes,"/home/pradeep/RNA_HOME/fusarium/ballgown/rename/Salt_vs_NoSalt_gene_results_sig.tsv",sep="\t")












gene_matrix = texpr(bg)
write.table(gene_matrix,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_gexpr.tsv",sep="\t")

trans = structure(bg)$trans
trans = trans[order(start(trans))]
exons = structure(bg)$exon

write.table(structure(bg)$trans,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_trans.tsv",sep="\t")
write.table(structure(bg)$exon,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_exons.tsv",sep="\t")
write.table(trans,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_trans_ordered.tsv",sep="\t")
write.table(exons,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_exons_ordered.tsv",sep="\t")







bg_table = texpr(bg, 'all')
knowngenetranscripts <- GRanges(seqnames=trans$"seqnames",ranges=IRanges(start=trans$"start",end=trans$"end"),strand=trans$"strand",spliceR.isoform_id = trans$"id",spliceR.gene_id=bg_table[match(trans$"group", bg_table$"t_id"),"gene_id"])
head(knowngenetranscripts)
knownGeneExons <- GRanges(seqnames=exons$"seqnames",ranges=IRanges(start=exons$"start",end=exons$"end"),strand=exons$"strand")
,spliceR.isoform_id = exons$"id",spliceR.gene_id=knowngenetranscripts[match(knowngenetranscripts$"spliceR.isoform_id", exons$"id"),"spliceR.gene_id"]
 
knowngenetranscripts <- GRanges(seqnames=results_transcripts$"chr",ranges=IRanges(start=results_transcripts$"start",end=results_transcripts$"end"),strand=results_transcripts$"strand",spliceR.isoform_id = results_transcripts$"t_name",spliceR.gene_id=results_transcripts$"gene_id")

library(spliceR)
bg_table = texpr(bg, 'all')
bg_fpkm = texpr(bg, meas = "FPKM")
head(bg_table)
write.table(bg_table,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_bgtable.tsv",sep="\t")
write.table(bg_fpkm,"/home/pradeep/RNA_HOME/fusarium/alternative_splicing/Salt_vs_NoSalt_bgfpkm.tsv",sep="\t")


splicerlist <- SpliceRList(structure(bg),exons, assembly_id = "pp", source = "fusarium", conditions = c("Salt","NoSalt"))
conditions(splicerlist)
exons(splicerlist)
dim(splicerlist)
splicer <- spliceR(splicerlist, compareTo = "preTranscript", filters = c("expressedGenes","geneOK", "isoOK", "expressedIso", "isoClass"), useProgressBar=F )
splicer <- spliceR(splicerlist)

# Exit the R session
# quit(save="no")