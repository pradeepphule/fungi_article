find out predicted 39 effector annotation using uniprotkb
then check for heatmap()
##########################################################################################################
#############secreated proteases identified by blastp of all transcripts against database sequences##############
##########################################################################################################
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/secreated_protein_new/FunSec-master/output/FunSec_Output/Final" 
setwd(working_dir)
secreated_protein = read.delim("seq_id.txt")
head(secreated_protein)
nrow(secreated_protein)
head(trans_annotation)
secreated_protein_annotation <- merge(secreated_protein,trans_annotation,by.x=c("transcript_id"),by.y=c("SeqName"))
head(secreated_protein_annotation)
nrow(secreated_protein_annotation)
secreated_protein_annotation$GOTerms<- NULL
secreated_protein_annotation$Enzymes<- NULL
secreated_protein_annotation$Pathway<- NULL
write.table(secreated_protein_annotation,"secreated_protein_annotation.tsv",sep="\t")
##########################################################################################################
#############effector proteins identified by blastp of all transcripts against database sequences##############
##########################################################################################################
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/secreated_protein_new/FunSec-master/output/FunSec_Output/Final/online_effectorP_result" 
setwd(working_dir)
effector_table = read.delim("predicted39_effector_id.txt")
nrow(effector_table)

#making a column of all significantly logfc2 id by combining all four conditions
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
sig_tid
write.table(sig_tid,"significant_id.tsv",sep="\t")

#merging all transporter containing transcripts with significant id 
sig_tid_effector_table <- merge(effector_table, sig_tid,by.x=c("t_id"),by.y=c("sig_tid"))
head(sig_tid_effector_table)
nrow(sig_tid_effector_table)

sig_tid_effector_table_tid<- unique(sig_tid_effector_table$t_id)

write.table(sig_tid_effector_table_tid,"significant_effector_table_supplementary.tsv",sep="\t")

bg_sig_tid_effector = subset(bg_all_filt, "t_name %in% sig_tid_effector_table_tid")
bg_sig_tid_effector

effector_transcript <- as.data.frame(texpr(bg_sig_tid_effector, meas='FPKM'))
effector_transcript


indices_effector <- match( rownames(effector_transcript), texpr(bg_sig_tid_effector, 'all')$t_id)
effector_names <- texpr(bg_sig_tid_effector, 'all')$t_name[indices_effector]
rownames(effector_transcript) <- effector_names
colnames(effector_transcript) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
effector_transcript
data_columns=c(1:12)
min_nonzero=1
effector_transcript
library(ComplexHeatmap)
library(circlize)

write.table(effector_transcript,"edit_effector_expression_transcript.tsv",sep="\t")


effector_transcript = log2(effector_transcript[,data_columns]+min_nonzero)

#transform data by row/gene centering
effector_transcript <- effector_transcript - matrix(rep(apply(effector_transcript,1,mean),ncol(effector_transcript)),ncol=ncol(effector_transcript))
head(effector_transcript)
#What are the minimum and maximum FPKM values for a particular library?
min(effector_transcript[,])
max(effector_transcript[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))


#sig_tid_cazy_table<-sig_tid_cazy_table[match(rownames(cazy_transcript),sig_tid_cazy_table$t_id),]


Heatmap(as.matrix(effector_transcript), name = "Expression",km = 3, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        ,top_annotation = ha, top_annotation_height = unit(1, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,show_row_names = TRUE,
        column_title = "Expression of Candidate effector encoding transcripts", 
        column_title_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 8),
        column_names_max_height = unit(2,"cm"),width = unit(5, "cm"))
####################################################
##########################################################################################################
#############cazy proteins identified by blastp of all transcripts against database sequences##############
##########################################################################################################
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/secreated_protein_new/FunSec-master/output/FunSec_Output/Final/dbCAN_output_only_secreated_protein" 
setwd(working_dir)
cazy_table = read.delim("input_cazy_cassified.txt")
nrow(cazy_table)
#cazy_table = unique(cazy_table,by=c("t_id"))
#nrow(cazy_table)
#making a column of all significantly logfc2 id by combining all four conditions
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
sig_tid
write.table(sig_tid,"significant_id.tsv",sep="\t")

#merging all transporter containing transcripts with significant id 
sig_tid_cazy_table <- merge(cazy_table, sig_tid,by.x=c("t_id"),by.y=c("sig_tid"))
head(sig_tid_cazy_table)
nrow(sig_tid_cazy_table)

sig_tid_cazy_table_tid<- unique(sig_tid_cazy_table$t_id)

sig_tid_cazy_table$t_id<-make.names(sig_tid_cazy_table$t_id, unique = TRUE)

write.table(sig_tid_cazy_table,"significant_cazy_table_supplementary.tsv",sep="\t")



bg_sig_tid_cazy = subset(bg_all_filt, "t_name %in% sig_tid_cazy_table_tid")
bg_sig_tid_cazy

cazy_transcript <- as.data.frame(texpr(bg_sig_tid_cazy, meas='FPKM'))
cazy_transcript


indices_cazy <- match( rownames(cazy_transcript), texpr(bg_sig_tid_cazy, 'all')$t_id)
cazy_names <- texpr(bg_sig_tid_cazy, 'all')$t_name[indices_cazy]
rownames(cazy_transcript) <- cazy_names
colnames(cazy_transcript) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
cazy_transcript
data_columns=c(1:12)
min_nonzero=1
cazy_transcript
library(ComplexHeatmap)
library(circlize)

write.table(cazy_transcript,"edit_cazy_expression_transcript.tsv",sep="\t")
cazy_transcript = read.table("edit_cazy_expression_transcript_edited.tsv",sep="\t",header  = TRUE)
#row.names(cazy_transcript) = cazy_transcript$t_id
rownames(cazy_transcript) <- make.names(cazy_transcript$X, unique = TRUE)
cazy_transcript$X<-NULL
cazy_transcript
cazy_transcript = log2(cazy_transcript[,data_columns]+min_nonzero)




# transform data by row/gene centering
cazy_transcript <- cazy_transcript - matrix(rep(apply(cazy_transcript,1,mean),ncol(cazy_transcript)),ncol=ncol(cazy_transcript))
head(cazy_transcript)
#What are the minimum and maximum FPKM values for a particular library?
min(cazy_transcript[,])
max(cazy_transcript[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))


sig_tid_cazy_table<-sig_tid_cazy_table[match(rownames(cazy_transcript),sig_tid_cazy_table$t_id),]


Heatmap(as.matrix(cazy_transcript), name = "Expression",km = 3, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        ,top_annotation = ha, top_annotation_height = unit(1, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,show_row_names = TRUE,
        column_title = "Expression of Cazyme encoding transcripts", 
        column_title_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 8),
        column_names_max_height = unit(2,"cm"),width = unit(5, "cm"))+Heatmap(as.matrix(sig_tid_cazy_table$categories),
                                                                              name="Categories",col = c("steelblue", "sienna2","seagreen3","yellowgreen"), 
                                                                              column_names_gp = gpar(fontsize = 8),
                                                                              column_names_max_height = unit(2,"cm"),
                                                                              width = unit(7, "mm"))
##########################################################################################################
#############transporter proteins identified by blastp of all transcripts against database sequences##############
##########################################################################################################
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/Fungal_cytochrome_P450 database" 
setwd(working_dir)
cytochrome_table = read.delim("cytochromeP450_result.txt")
nrow(cytochrome_table)
cytochrome_table = unique(cytochrome_table,by=c("t_id"))
nrow(cytochrome_table)
#making a column of all significantly logfc2 id by combining all four conditions
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
sig_tid
write.table(sig_tid,"significant_id.tsv",sep="\t")
#merging all transporter containing transcripts with significant id 
sig_tid_cytochrome_table <- merge(cytochrome_table, sig_tid,by.x=c("t_id"),by.y=c("sig_tid"))
head(sig_tid_cytochrome_table)
nrow(sig_tid_cytochrome_table)

write.table(sig_tid_cytochrome_table,"significant_cytochrome_table_supplementary.tsv",sep="\t")

sig_tid_cytochrome_table_tid<- sig_tid_cytochrome_table$t_id

bg_sig_tid_cytochrome = subset(bg_all_filt, "t_name %in% sig_tid_cytochrome_table_tid")
bg_sig_tid_cytochrome

cytochrome_transcript <- as.data.frame(texpr(bg_sig_tid_cytochrome, meas='FPKM'))
cytochrome_transcript
indices_cytochrome <- match( rownames(cytochrome_transcript), texpr(bg_sig_tid_cytochrome, 'all')$t_id)
cytochrome_names <- texpr(bg_sig_tid_cytochrome, 'all')$t_name[indices_cytochrome]
rownames(cytochrome_transcript) <- cytochrome_names
colnames(cytochrome_transcript) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
data_columns=c(1:12)
min_nonzero=1
cytochrome_transcript
library(ComplexHeatmap)
library(circlize)
cytochrome_transcript = log2(cytochrome_transcript[,data_columns]+min_nonzero)
head(cytochrome_transcript)
# transform data by row/gene centering
cytochrome_transcript <- cytochrome_transcript - matrix(rep(apply(cytochrome_transcript,1,mean),ncol(cytochrome_transcript)),ncol=ncol(cytochrome_transcript))
head(cytochrome_transcript)
#What are the minimum and maximum FPKM values for a particular library?
min(cytochrome_transcript[,])
max(cytochrome_transcript[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))
Heatmap(as.matrix(cytochrome_transcript), name = "Expression",km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        ,top_annotation = ha, top_annotation_height = unit(1, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,show_row_names = TRUE,row_names_gp = gpar(fontsize = 5),
        column_title = "Expression of CytP450 encoding transcripts", 
        column_title_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 8),
        column_names_max_height = unit(2,"cm"),width = unit(5, "cm"))







##########################################################################################################
#############transporter proteins identified by blastp of all transcripts against database sequences##############
##########################################################################################################
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/TCDB" 
setwd(working_dir)
transporter_blast_table = read.delim("final_supplementary.txt")
transporter_table = read.delim("final_supplementary_change.txt")
nrow(transporter_table)#3453
head(transporter_table)
transporter_table_unique = unique(transporter_table,by=c("t_id"))
nrow(transporter_table_unique)#2726
head(transporter_table_unique)
#unique_tid<-as.matrix(unique(transporter_table$t_id))
#colnames(unique_tid)<-"unique_tid"
#head(unique_tid)
#nrow(unique_tid)
#transporter_table = merge(unique_tid,transporter_table,by.x=c("unique_tid"),by.y=c("t_id"))
#head(transporter_table)
#nrow(transporter_table)
#making a column of all significantly logfc2 id by combining all four conditions
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
nrow(sig_tid)#2698
head(sig_tid)
#merging all transporter containing transcripts with significant id 
sig_tid_transporter_table <- merge(transporter_table_unique, sig_tid,by.x=c("t_id"),by.y=c("sig_tid"))
head(sig_tid_transporter_table)
nrow(sig_tid_transporter_table)#501

# mapped information from stringtie is ordered according to core geneid and transcriptid of stringtie
#refmap_stringtie = read.delim("merged.stringtie_merged.gtf.refmap")
#colnames(refmap_stringtie)
#head(refmap_stringtie)

#refmap_stringtie_geneid = read.delim("stringtie_geneid_transcriptid.refmap")
#colnames(refmap_stringtie_geneid)
#head(refmap_stringtie_geneid)
#trans_transporter_id <- merge(transporter_table,refmap_stringtie_geneid,by.x=c("t_id"),by.y=c("t_name"),all.x=TRUE)

trans_transporter_id <- merge(transporter_blast_table, tmap_stringtie,by.x=c("t_id"),by.y=c("transcript_id"),all.x=TRUE)
head(trans_transporter_id)
nrow(trans_transporter_id)
sig_tid_transporter_table_match_id <- merge(sig_tid_transporter_table, tmap_stringtie,by.x=c("t_id"),by.y=c("transcript_id"),all.x=TRUE)
head(sig_tid_transporter_table_match_id)
nrow(sig_tid_transporter_table_match_id)
write.table(trans_transporter_id,"trnascripts_to_transporter_id_match_table.tsv",sep="\t")
write.table(sig_tid_transporter_table_match_id,"sig_tid_trnascripts_to_transporter_id_match_table.tsv",sep="\t")



# split the 501 transcripts in four condition
#head(sig1_logfc2)
#nrow(sig1_logfc2)
#sig1_logfc2_transporter_table <- match(sig1_logfc2$feature, sig_tid_transporter_table$t_id)

sig1_logfc2_transporter_table <- merge(sig1_logfc2, sig_tid_transporter_table,by.x=c("feature"),by.y=c("t_id"))
head(sig1_logfc2_transporter_table)
nrow(sig1_logfc2_transporter_table)#126
#sig1_logfc2_transporter_table <- merge(sig1_logfc2_transporter_table, tmap_stringtie,by.x=c("feature"),by.y=c("transcript_id"),all.x=TRUE)
write.table(sig1_logfc2_transporter_table,"sig1_logfc2_transporter_table.tsv",sep="\t")


sig2_logfc2_transporter_table <- merge(sig2_logfc2, sig_tid_transporter_table,by.x=c("feature"),by.y=c("t_id"))
head(sig2_logfc2_transporter_table)
nrow(sig2_logfc2_transporter_table)#245
#sig2_logfc2_transporter_table <- merge(sig2_logfc2_transporter_table, tmap_stringtie,by.x=c("feature"),by.y=c("transcript_id"),all.x=TRUE)
write.table(sig2_logfc2_transporter_table,"sig2_logfc2_transporter_table.tsv",sep="\t")

sig3_logfc2_transporter_table <- merge(sig3_logfc2, sig_tid_transporter_table,by.x=c("feature"),by.y=c("t_id"))
head(sig3_logfc2_transporter_table)
nrow(sig3_logfc2_transporter_table)#280
#sig3_logfc2_transporter_table <- merge(sig3_logfc2_transporter_table, tmap_stringtie,by.x=c("feature"),by.y=c("transcript_id"),all.x=TRUE)
write.table(sig3_logfc2_transporter_table,"sig3_logfc2_transporter_table.tsv",sep="\t")

sig4_logfc2_transporter_table <- merge(sig4_logfc2, sig_tid_transporter_table,by.x=c("feature"),by.y=c("t_id"))
head(sig4_logfc2_transporter_table)
nrow(sig4_logfc2_transporter_table)#75
#sig4_logfc2_transporter_table <- merge(sig4_logfc2_transporter_table, tmap_stringtie,by.x=c("feature"),by.y=c("transcript_id"),all.x=TRUE)
write.table(sig4_logfc2_transporter_table,"sig4_logfc2_transporter_table.tsv",sep="\t")



sig_tid_transporter_table = arrange(sig_tid_transporter_table,Transporter_class)
head(sig_tid_transporter_table)
nrow(sig_tid_transporter_table)#501
write.table(sig_tid_transporter_table,"significant_final_supplementary.tsv",sep="\t")

Transporter_class<- sig_tid_transporter_table$Transporter_class

Transporter_class_tid<- sig_tid_transporter_table$t_id

bg_sig_tid_transporter = subset(bg_all_filt, "t_name %in% Transporter_class_tid")
bg_sig_tid_transporter#ballgown instance with 501 transcripts and 12 samples

transporter_transcript <- as.data.frame(texpr(bg_sig_tid_transporter, meas='FPKM'))
transporter_transcript
indices_transporter <- match( rownames(transporter_transcript), texpr(bg_sig_tid_transporter, 'all')$t_id)
transporter_names <- texpr(bg_sig_tid_transporter, 'all')$t_name[indices_transporter]
rownames(transporter_transcript) <- transporter_names
colnames(transporter_transcript) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
data_columns=c(1:12)
min_nonzero=1
nrow(transporter_transcript)#501
head(transporter_transcript)
library(ComplexHeatmap)
library(circlize)


transporter_transcript = log2(transporter_transcript[,data_columns]+min_nonzero)
head(transporter_transcript)
# transform data by row/gene centering
transporter_transcript <- transporter_transcript - matrix(rep(apply(transporter_transcript,1,mean),ncol(transporter_transcript)),ncol=ncol(transporter_transcript))
head(transporter_transcript)
#What are the minimum and maximum FPKM values for a particular library?
min(transporter_transcript[,])
max(transporter_transcript[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))

transporter_table_significant_logfc2<-sig_tid_transporter_table[match(rownames(transporter_transcript),sig_tid_transporter_table$t_id),]
nrow(transporter_table_significant_logfc2)#501
head(transporter_table_significant_logfc2)
write.table(transporter_table_significant_logfc2,"transporter_table_significant_logfc2.tsv",sep="\t")

Heatmap(as.matrix(transporter_transcript), name = "Expression",km = 5,
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        ,top_annotation = ha, top_annotation_height = unit(1, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,show_row_names = TRUE,
       row_names_gp = gpar(fontsize = 6),
        column_title = "Expression of transporter encoding transcripts", 
        column_title_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 8),
        column_names_max_height = unit(2,"cm"),width = unit(5, "cm"))+Heatmap(as.matrix(transporter_table_significant_logfc2$Transporter_class),
                                        name="Transporter_class",col = c("blue", "thistle4", "sienna", "Orchid","sienna2","seagreen3","yellowgreen"), 
                                        column_names_gp = gpar(fontsize = 6),
                                        column_names_max_height = unit(2,"cm"),
                                        width = unit(7, "mm"))


##########################################################################################################
#Clean up workspace - i.e. delete variable created by the graphics demo
#rm(list = ls(all = TRUE))
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/transction_factors" 
setwd(working_dir)
trans_interproid = read.delim("blast2go_export_interproscan_id.txt")
fungal_interproid = read.delim("fungal_transcrition_factors_interproscanid.txt")
trans_topinteroid = read.delim("trans_interproscan_topid.txt")

DF3 <- merge(trans_topinteroid, fungal_interproid,by.x=c("First_InterPro_Term"),by.y=c("InterPro_Term"))
head(DF3)
write.table(DF3,"trnascripts_to_t.fators_match_table.tsv",sep="\t")

sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
sig_tid

DF3_sig_tid <- merge(sig_tid,DF3 ,by.x=c("sig_tid"),by.y=c("transcript_name"))
DF3_sig_tid
trans_factor_tid<-DF3_sig_tid$sig_tid
bg_trans_factor = subset(bg_all_filt, "t_name %in% trans_factor_tid")
bg_trans_factor
trans_factor <- as.data.frame(texpr(bg_trans_factor, meas='FPKM'))
trans_factor

indices_trans_factor <- match( rownames(trans_factor), texpr(bg_all_filt, 'all')$t_id)
trans_factor_names <- texpr(bg_all_filt, 'all')$t_name[indices_trans_factor]
rownames(trans_factor) <- trans_factor_names

colnames(trans_factor) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
data_columns=c(1:12)
min_nonzero=1
trans_factor
library(ComplexHeatmap)
library(circlize)
trans_factor = log2(trans_factor[,data_columns]+min_nonzero)
head(trans_factor)
# transform data by row/gene centering
trans_factor <- trans_factor - matrix(rep(apply(trans_factor,1,mean),ncol(trans_factor)),ncol=ncol(trans_factor))
head(trans_factor)
#What are the minimum and maximum FPKM values for a particular library?
min(trans_factor[,])
max(trans_factor[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))

DF3_sig_tid_significant_logfc2<-DF3_sig_tid[match(rownames(trans_factor),DF3_sig_tid$sig_tid),]

write.table(DF3_sig_tid_significant_logfc2,"DF3_sig_tid_significant_logfc2_table.tsv",sep="\t")


Heatmap(as.matrix(trans_factor), name = "Expression",km = 5,
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        ,top_annotation = ha, top_annotation_height = unit(1, "mm")
        ,cluster_columns = FALSE, cluster_rows = TRUE
        ,column_title = "Expression of transcription factors encoding transcripts"
        ,row_dend_width = unit(10, "mm"), 
        row_names_max_width = unit(10, "cm"),
        row_names_gp = gpar(fontsize = 5),
        column_title_gp = gpar(fontsize = 12),
        column_names_max_height = unit(2,"cm"),width = unit(10, "cm"))+Heatmap(as.matrix(DF3_sig_tid_significant_logfc2$TF_Family_Name), 
                                                                               name="TF_Family",
                                                                               col = c("steelblue","springgreen" ,"cyan","red1","yellowgreen","yellow" ,"darkblue","chocolate","violetred" ,"darkkhaki","darkgrey"),
                                                                               column_names_gp = gpar(fontsize = 9),column_names_max_height = unit(2,"cm"),width = unit(10, "mm"))


##########################################################################################################
#############find out secondary metabolite mapped genes with stringtie genes and transcripts##############
##########################################################################################################
#rm(list = ls(all = TRUE)) first execute code till bg_all_filt
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/antismash_smurf" 
setwd(working_dir)
# antismash core geneid , type of SM and cluster id information is extracted
core_sm_id = read.delim("smid_clusterid_smtype.file")
colnames(core_sm_id)
head(core_sm_id)
# mapped information from stringtie is ordered according to core geneid and transcriptid of stringtie
refmap_stringtie = read.delim("merged.stringtie_merged.gtf.refmap")
colnames(refmap_stringtie)
refmap_stringtie$class_code<- NULL
refmap_stringtie$ref_id<- NULL
head(refmap_stringtie)
trans_SM_id <- merge(core_sm_id, refmap_stringtie,by.x=c("ref_gene_id"),by.y=c("ref_gene_id"))
head(trans_SM_id)
write.table(trans_SM_id,"trnascripts_to_sm_coregeneid_match_table.tsv",sep="\t")
#order match table according to refgeneid-antismashclusterid-type-stringtiegeneid-transcript manually
mapped_sm_id = read.delim("trnascripts_to_sm_coregeneid_match_table_ordered.txt")
head(mapped_sm_id)
#run the code from transcript_assign_main.r upto bg_all table
# retrive transcript level fpkm values of transcript ids.from ballgown subset and make matrix of these ids.
#t_n significant by p value 0.05
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
head(sig_tid)

sig_tid_mapped_sm_id <- merge(mapped_sm_id, sig_tid,by.x=c("ST_transcript_id"),by.y=c("sig_tid"))
sig_tid_mapped_sm_id
write.table(sig_tid_mapped_sm_id,"mapped_sm_id_significant_logfc2_table.tsv",sep="\t")
mapped_sm_id = arrange(mapped_sm_id,Type)
mapped_sm_id
sm_stringtieid<-mapped_sm_id$ST_transcript_id
unique(sm_stringtieid)
bg_sm_transcripts = subset(bg_all_filt, "t_name %in% sm_stringtieid")
bg_sm_transcripts
sm_transcripts <- as.data.frame(texpr(bg_sm_transcripts, meas='FPKM'))
sm_transcripts
indices <- match( rownames(sm_transcripts), texpr(bg_all_filt, 'all')$t_id)
sm_names <- texpr(bg_all_filt, 'all')$t_name[indices]
rownames(sm_transcripts) <- sm_names
colnames(sm_transcripts) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
sm_transcripts
data_columns=c(1:12)
min_nonzero=1
sm_transcripts
library(ComplexHeatmap)
library(circlize)
sm_transcripts = log2(sm_transcripts[,data_columns]+min_nonzero)
head(sm_transcripts)
# transform data by row/gene centering
sm_transcripts <- sm_transcripts - matrix(rep(apply(sm_transcripts,1,mean),ncol(sm_transcripts)),ncol=ncol(sm_transcripts))
head(sm_transcripts)
#What are the minimum and maximum FPKM values for a particular library?
min(sm_transcripts[,])
max(sm_transcripts[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))


unique(mapped_sm_id$ST_transcript_id)
mapped_sm_id_match_type<-mapped_sm_id[match(rownames(sm_transcripts),mapped_sm_id$ST_transcript_id),]

write.table(mapped_sm_id_match_type,"mapped_sm_id_match_type.tsv",sep="\t")

Heatmap(as.matrix(sm_transcripts), name = "Expression"
        ,km = 3, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(3, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,
        column_title = "Expression of core secondary Metabolite transcripts",
        column_title_gp = gpar(fontsize = 12),column_names_max_height = unit(2,"cm"),
        width = unit(15, "cm"))+Heatmap(as.matrix(mapped_sm_id_match_type$Type), 
                                        name="SM_Type",
                                        col = c("steelblue", "Orchid","sienna2","seagreen3","yellowgreen","chocolate"),
                                        column_names_max_height = unit(2,"cm"),width = unit(10, "mm"))

##########################################################################################################
##########
#rm(list = ls(all = TRUE)) first execute code till bg_all_filt
working_dir = "/home/pradeep/RNA_HOME/fusarium/stringtie/de_novo/antismash_smurf" 
setwd(working_dir)
# antismash core geneid , type of SM and cluster id information is extracted
core_sm_id = read.delim("smid_clusterid_smtype.file")
colnames(core_sm_id)
head(core_sm_id)
# mapped information from stringtie is ordered according to core geneid and transcriptid of stringtie
tmap_stringtie = read.delim("mapped_file.tmap")
head(tmap_stringtie)
trans_SM_id <- merge(core_sm_id, tmap_stringtie,by.x=c("ref_gene_id"),by.y=c("ref_gene_id"))
head(trans_SM_id)
#run the code from transcript_assign_main.r upto bg_all table
# retrive transcript level fpkm values of transcript ids.from ballgown subset and make matrix of these ids.
#t_n significant by p value 0.05
sig_tid<-data.frame(t_n_logfc2)#from main script combination of all significant transcripts
colnames(sig_tid)<-"sig_tid"
head(sig_tid)
sig_tid_trans_SM_id <- merge(trans_SM_id, sig_tid,by.x=c("transcript_id"),by.y=c("sig_tid"))
sig_tid_trans_SM_id
write.table(sig_tid_trans_SM_id,"tmap_based_sm_id_significant_logfc2_table.tsv",sep="\t")
sig_tid_trans_SM_id = arrange(sig_tid_trans_SM_id,Type)
sig_tid_trans_SM_id

sm_stringtieid<-sig_tid_trans_SM_id$transcript_id
unique(sm_stringtieid)
bg_sm_transcripts = subset(bg_all_filt, "t_name %in% sm_stringtieid")
bg_sm_transcripts

sm_transcripts <- as.data.frame(texpr(bg_sm_transcripts, meas='FPKM'))
sm_transcripts
indices <- match( rownames(sm_transcripts), texpr(bg_all_filt, 'all')$t_id)
sm_names <- texpr(bg_all_filt, 'all')$t_name[indices]
rownames(sm_transcripts) <- sm_names
colnames(sm_transcripts) <- c("Salt.r1","Salt.r2","Salt.r3","NoSalt.r1","NoSalt.r2","NoSalt.r3","Bac.r1","Bac.r2","Bac.r3","Metabac.r1","Metabac.r2","Metabac.r3")
sm_transcripts
data_columns=c(1:12)
min_nonzero=1
sm_transcripts
library(ComplexHeatmap)
library(circlize)
sm_transcripts = log2(sm_transcripts[,data_columns]+min_nonzero)
head(sm_transcripts)
# transform data by row/gene centering
sm_transcripts <- sm_transcripts - matrix(rep(apply(sm_transcripts,1,mean),ncol(sm_transcripts)),ncol=ncol(sm_transcripts))
head(sm_transcripts)
#What are the minimum and maximum FPKM values for a particular library?
min(sm_transcripts[,])
max(sm_transcripts[,])
type <- c("r1","r2","r3","r1","r2","r3","r1","r2","r3","r1","r2","r3")
ha = HeatmapAnnotation(df = data.frame(type = type))




unique(sig_tid_trans_SM_id$transcript_id)
sig_tid_trans_SM_id_match_type<-sig_tid_trans_SM_id[match(rownames(sm_transcripts),sig_tid_trans_SM_id$transcript_id),]

#write.table(sig_tid_trans_SM_id_match_type,"sig_tid_trans_SM_id_match_type.tsv",sep="\t")

Heatmap(as.matrix(sm_transcripts), name = "Expression"
        ,km = 3, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(3, "mm"), 
        cluster_columns = FALSE, cluster_rows = TRUE,
        column_title = "Expression of core secondary metabolite transcripts",
        column_title_gp = gpar(fontsize = 12),column_names_max_height = unit(2,"cm"),
        width = unit(15, "cm"))+Heatmap(as.matrix(sig_tid_trans_SM_id_match_type$type), 
                                        name="SM_Type",
                                        col = c("steelblue", "Orchid","sienna2","seagreen3","yellowgreen","chocolate"),
                                        column_names_max_height = unit(2,"cm"),width = unit(10, "mm"))
