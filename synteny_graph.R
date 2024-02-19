
library(genoPlotR)
df2 <- data.frame(name=c("ctg27_orf42", "ctg27_orf43", "ctg27_orf44",
                         "ctg27_orf45", "ctg27_orf46", "ctg27_orf47", 
                         "ctg27_orf48", "ctg27_orf49", "ctg27_orf50",
                         "ctg27_orf51")
                  ,start=c(11159,13005,17586,
                           20001,32436,36093,
                           38221,40363,44884,
                           49642)
                  ,end=c(11439,14504,18513,
                         32068,35791,37531,
                         38966,44708,48041,
                         50412)
                  ,strand=c(1,1,1,
                            -1,1,-1,
                            1,1,1,
                            1)
                  ,col=c("blue", "grey", "red",
                         "blue", "grey", "red",
                         "blue", "grey", "red",
                         "blue" ))

dna_seg2 <- dna_seg(df2)


df1 <- data.frame(name=c("KIL85384", "KIL85385", "KIL85386", "KIL85387", "KIL85388",
                         "KIL85389", "KIL85390", "KIL85391", "KIL85392", "KIL85393",
                         "KIL85394", "KIL85395", "KIL85396", "KIL85397", "KIL85398", 
                         "KIL85399", "KIL85400", "KIL85401", "KIL85402", "KIL85403",
                         "KIL85404", "KIL85405", "KIL85406", "KIL85407", "KIL85408", "KIL85409")
                  ,start=c(291146,293919,296024,298754,302767,
                           305097,309148,311940,314364,319725,
                           332528,334759,335813,337940,340082,
                           342197,350968,353224,356289,358337,
                           360083,361352,363492,366434,366928,369090)
                  ,end=c(293091,295664,297661,300445,304278,
                         307519,311485,313154,319258,331629,
                         333787,335512,337150,338641,341812,
                         350728,352158,355007,357017,359495,
                         361006,362893,365112,366783,368127,370726
                  )
                  ,strand=c(1,1,-1,1,1,
                            -1,1,1,1,-1,
                            1,1,-1,1,1,
                            1,-1,-1,-1,-1,
                            1,1,-1,1,-1,-1)
                  ,col=c("blue", "grey", "red","blue", "grey", 
                         "red","blue", "grey", "red","blue",
                         "grey", "red","blue", "grey", "red",
                         "blue", "grey", "red","blue", "grey",
                         "red","blue", "grey", "red","blue", "grey"))

dna_seg1 <- dna_seg(df1)

df3 <- data.frame(name=c("EKJ71904","EKJ71905","EKJ71906","EKJ71907","EKJ71908",
                         "EKJ71909","EKJ71910","EKJ71911","EKJ71912","EKJ71913",
                         "EKJ71914","EKJ71915","EKJ71916","EKJ71917","EKJ71918",
                         "EKJ71919","EKJ71920","EKJ71921","EKJ71922","EKJ71923",
                         "EKJ71924","EKJ71925")
                  ,start=c(39393,42578,45538,48792,50365,
                           51849,54584,56672,69273,71345,
                           72285,74372,76333,78284,80762,
                           82810,88958,89569,92675,94582,
                           96185,103556)
                  ,end=c(40037,43681,46593,49954,51291,
                         53978,55894,68534,70532,72066,
                         73598,75073,78009,80088,82417,
                         83925,89441,91643,93964,95661,
                         97663,104045)
                  ,strand=c(1,-1,1,1,-1,
                            -1,-1,-1,1,1,
                            -1,1,1,1,1,
                            1,1,-1,1,1,
                            -1,1)
                  ,col=c("blue", "grey", "red","blue", "grey", 
                         "red","blue", "grey", "red","blue",
                         "grey", "red","blue", "grey", "red",
                         "blue", "grey", "red","blue", "grey",
                         "red","blue"))

dna_seg3 <- dna_seg(df3)

dna_segs <- list(dna_seg1, dna_seg2, dna_seg3)
dna_segs1 <- list(dna_seg1, dna_seg2)
df4 <- data.frame(start1=c(309148,311940,319725,332528,334759,335813,337940,340082,342197,342197,342197),
                  end1=c(311485,313154,331629,333787,335512,337150,338641,341812,350728,350728,350728),
  start2=c(11159,13005,20001,32436,32436,36093,38221,40363,40363,44884,49642),
                  end2=c(11439,14504,32068,35791,35791,37531,38966,44708,44708,48041,50412)
                  
           #       ,direction=c()
)
#in comparision object put two colors dark gery for differnt identities from clusterblast output. or you can also make some grades of colors
comparison1 <- comparison(df4)
df5 <- data.frame(start1=c(20001,32436,32436,36093,38221,40363,40363,44884,44884),
                  end1=c(32068,35791,35791,37531,38966,44708,44708,48041,48041),
                  start2=c(56672,69273,71345,72285,74372,76333,78284,80762,82810),
                  end2=c(68534,70532,72066,73598,75073,78009,80088,82417,83925)
          #        ,direction=c()
                  )                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
comparison2 <- comparison(df5)

names <- c("Huey", "Dewey", "Louie")
names(dna_segs) <- names

comparisons <- list(comparison1,comparison2)
comparisons1 <- list(comparison1)
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, main="Comparison of Huey, Dewey and Louie")
plot_gene_map(dna_segs=dna_segs1, comparisons=comparisons1, main="Comparison of Huey, Dewey and Louie")


