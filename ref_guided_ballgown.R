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
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/ref_guided/ballgown" 
setwd(working_dir)
pheno_data = read.csv("all_path.csv")



########################################################################################################################
bg_all = ballgown(samples=as.vector(pheno_data$path),pData=pheno_data, meas='all')
bg_all # ballgown instance with 25208 transcripts and 12 samples
bg_all_table = texpr(bg_all, 'all')
head(bg_all_table)
bg_all_filt = subset(bg_all,"rowVars(texpr(bg_all)) > 1", genomesubset=TRUE)#Low abundance genes were filtered by removing transcripts with a variance across the sample less than one. 
bg_all_filt # ballgown instance with 14234 transcripts and 12 samples