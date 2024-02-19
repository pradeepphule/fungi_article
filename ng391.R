library(genoPlotR)

#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))

df1 <- data.frame(name=c("no_tag_f_4659","no_tag_f_4660","no_tag_f_4661",
                         "no_tag_f_4662","no_tag_f_4663","no_tag_f_4664",
                         "hybrid PKS-NRPS","alpha/beta hydrolase","eEF-1B gamma subunit",
                         "no_tag_f_4668","no_tag_f_4669","MFS protein",
                         "no_tag_f_4671","no_tag_f_4672"),
                  start=c(144797,147488,153157,
                          157246,160142,164723,
                          167138,179948,182169,
                          183230,185358,187500,
                          189615,197906),
                  end=c(147236,149234,155329,
                        159258,161357,166168,
                        179049,181208,182929,
                        184568,186060,189231,
                        197550,198980),
                  strand=c(-1,1,-1,
                           1,1,1,
                           -1,1,1,
                           -1,1,1,
                           1,-1),
                  col=c("black", "black","black",
                        "black", "black","black",
                        "red", "green", "blue",
                        "black", "black", "yellow",
                        "black", "black"))
dna_seg1 <- dna_seg(df1)

df2 <- data.frame(name=c("ACS68551.1","ACS68552.1","ACS68553.1",
                         "ACS68554.1","ACS68555.1","ACS68556.1"),
                  start=c(1,1264,2466,
                          4157,16691,17452),
                  end=c(694,2002,3746,
                        16022,17191,19101),
                  strand=c(1,-1,-1,
                           1,-1,1),
                  col=c("black", "blue", "green",
                        "red", "black", "yellow"))
dna_seg2 <- dna_seg(df2)
dna_segs <- list(dna_seg1, dna_seg2)

names <- c("F. pp1-1 
NG-391 like gene cluster", "NG-391  
biosynthetic gene cluster")
names(dna_segs) <- names

df3 <- data.frame(start1=c(167138,179948,182169,187500),
                  end1=c(179049,181208,182929,189231),
                  start2=c(4157, 2466,1264,17452),
                  end2=c(16022, 3746,2002,19101)
                 # col=c("#67000D", "#08306B","#67000D", "#08306B")
                  )
comparison1 <- comparison(df3)
comparisons <- list(comparison1)
mid_pos <- middle(dna_segs[[1]])
annot <- annotation(x1=c(mid_pos[7],mid_pos[8],mid_pos[9],mid_pos[12]),
                    x2=c(NA,NA,NA,NA),
                    text=c(dna_segs[[1]]$name[7],dna_segs[[1]]$name[8],
                           dna_segs[[1]]$name[9],dna_segs[[1]]$name[12]),
                    rot=c(30,40,30,20), col=c("red","green","blue","yellow"))
                    
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
              annotations=annot, annotation_height=4.5,
             # main="Comparison of Huey, Dewey and Louie"
              )


