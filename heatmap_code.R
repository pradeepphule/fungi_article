library(ComplexHeatmap)
library(circlize)
#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))
sm_type <- read.csv("~/bin/Final_clcbioandabyss/analysis/sm_type", sep="")
  View(sm_type)
  sm_type1 <- read.csv("~/bin/Final_clcbioandabyss/analysis/sm_type1", sep="")
  View(sm_type1)
  blast_sm <- read.csv("~/bin/Final_clcbioandabyss/analysis/blast_sm", sep="")
  View(blast_sm)
 ht=Heatmap(as.matrix(sm_type1), name="Type",col = c("lightblue", "Brown",
                                                     "Peru", "lavender","Orchid","grey","cyan","coral","Darkcyan","Pink"),
            column_names_max_height = unit(0.5,"cm"))
 + Heatmap(as.matrix(blast_sm), name = "Identity", 
           col = colorRamp2(c(0,0.4,1),c("white", "white","orange")), 
           cluster_columns = FALSE, cluster_rows = FALSE,
           column_title = "Occurrence of SM genes",
           column_title_gp = gpar(fontsize = 10),column_names_max_height = unit(2,"cm"))
draw(ht, heatmap_legend_side = "left")


ht=Heatmap(as.matrix(mapped_sm_id$Type), name="Type",col = c("blue", "Brown","Peru", "Pink","Orchid","grey","cyan","coral"),column_names_max_height = unit(0.5,"cm"))+ Heatmap(as.matrix(sm_transcripts), name = "Expression",km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")), cluster_columns = TRUE, cluster_rows = TRUE, column_title = "Transcript level expression of core secondary Metabolite genes", column_title_gp = gpar(fontsize = 15),column_names_max_height = unit(2,"cm"))
draw(ht, heatmap_legend_side = "right")

ht= Heatmap(as.matrix(sm_transcripts), name = "Expression",km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),  cluster_columns = FALSE, cluster_rows = TRUE, column_title = "Transcript level expression of core secondary Metabolite genes", column_title_gp = gpar(fontsize = 15),column_names_max_height = unit(2,"cm"))+Heatmap(as.matrix(mapped_sm_id$Type), name="Type",col = c("blue", "Brown","Peru", "Pink","Orchid","grey","cyan","coral"),column_names_max_height = unit(0.5,"cm"))
draw(ht, heatmap_legend_side = "right")








