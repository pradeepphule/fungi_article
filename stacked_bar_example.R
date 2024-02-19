#https://stackoverflow.com/questions/36945185/display-color-in-ggplot-2-with-scale-fill-manual-and-scale-fill-discrete-in-stac
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(reshape2)
library(sciplot)
library(data.table)
library(lattice)

weedseed3<-structure(list(Crop = structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                  2L, 2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
                                  4L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 
                                  3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
                                  4L, 4L, 4L, 4L, 4L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                  1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 
                                  3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
                                  4L, 4L, 4L, 4L, 4L), .Label = c("alfalfa", "corn", "oat", "soybean"
                                  ), class = "factor"), Rot = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 
                                                                          1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                          1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                          2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                          2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                                                          3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                                                          3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
                                                                          3L, 3L, 3L, 3L, 3L, 3L), .Label = c("2-year", "3-year", "4-year"
                                                                          ), class = "factor"), Rot.Herb = structure(c(3L, 3L, 3L, 3L, 
                                                                                                                       3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 13L, 13L, 13L, 13L, 13L, 13L, 
                                                                                                                       14L, 14L, 14L, 14L, 14L, 14L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 
                                                                                                                       6L, 6L, 6L, 6L, 9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L, 
                                                                                                                       10L, 15L, 15L, 15L, 15L, 15L, 15L, 16L, 16L, 16L, 16L, 16L, 16L, 
                                                                                                                       1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 7L, 7L, 7L, 7L, 
                                                                                                                       7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 11L, 11L, 11L, 11L, 11L, 11L, 
                                                                                                                       12L, 12L, 12L, 12L, 12L, 12L, 17L, 17L, 17L, 17L, 17L, 17L, 18L, 
                                                                                                                       18L, 18L, 18L, 18L, 18L), .Label = c("A4-conv", "A4-low", "C2-conv", 
                                                                                                                                                            "C2-low", "C3-conv", "C3-low", "C4-conv", "C4-low", "O3-conv", 
                                                                                                                                                            "O3-low", "O4-conv", "O4-low", "S2-conv", "S2-low", "S3-conv", 
                                                                                                                                                            "S3-low", "S4-conv", "S4-low"), class = "factor"), species = c("AMATA", 
                                                                                                                                                                                                                           "CHEAL", "ABUTH", "Others", "SETSP", "SOLPT", "CHEAL", "AMATA", 
                                                                                                                                                                                                                           "ABUTH", "Others", "SETSP", "SOLPT", "AMATA", "CHEAL", "SETSP", 
                                                                                                                                                                                                                           "SOLPT", "Others", "ABUTH", "CHEAL", "AMATA", "SETSP", "SOLPT", 
                                                                                                                                                                                                                           "Others", "ABUTH", "AMATA", "CHEAL", "SOLPT", "SETSP", "Others", 
                                                                                                                                                                                                                           "ABUTH", "CHEAL", "AMATA", "SETSP", "SOLPT", "Others", "ABUTH", 
                                                                                                                                                                                                                           "AMATA", "CHEAL", "SOLPT", "Others", "SETSP", "ABUTH", "CHEAL", 
                                                                                                                                                                                                                           "AMATA", "SETSP", "ABUTH", "Others", "SOLPT", "CHEAL", "AMATA", 
                                                                                                                                                                                                                           "Others", "ABUTH", "SETSP", "SOLPT", "CHEAL", "AMATA", "SETSP", 
                                                                                                                                                                                                                           "SOLPT", "ABUTH", "Others", "AMATA", "CHEAL", "SETSP", "Others", 
                                                                                                                                                                                                                           "SOLPT", "ABUTH", "AMATA", "CHEAL", "SETSP", "Others", "SOLPT", 
                                                                                                                                                                                                                           "ABUTH", "CHEAL", "SETSP", "AMATA", "Others", "ABUTH", "SOLPT", 
                                                                                                                                                                                                                           "SETSP", "AMATA", "CHEAL", "Others", "ABUTH", "SOLPT", "AMATA", 
                                                                                                                                                                                                                           "CHEAL", "SETSP", "SOLPT", "ABUTH", "Others", "AMATA", "CHEAL", 
                                                                                                                                                                                                                           "SETSP", "Others", "ABUTH", "SOLPT", "AMATA", "CHEAL", "SETSP", 
                                                                                                                                                                                                                           "SOLPT", "ABUTH", "Others", "CHEAL", "AMATA", "SOLPT", "SETSP", 
                                                                                                                                                                                                                           "Others", "ABUTH"), meandensity = c(782.0372807, 539.727819025, 
                                                                                                                                                                                                                                                               0, 0, 0, 0, 674.0338685, 638.871908975, 0, 0, 0, 0, 2993.61661595, 
                                                                                                                                                                                                                                                               1920.306528775, 36.0191411, 36.0191411, 35.224044875, 0, 2435.059165225, 
                                                                                                                                                                                                                                                               1829.35511415, 191.16853105, 78.12902155, 39.064510775, 0, 2717.26455395, 
                                                                                                                                                                                                                                                               1997.411972825, 818.117344125, 570.7718258, 106.5769257, 0, 2777.305886425, 
                                                                                                                                                                                                                                                               1687.689335975, 970.704059675, 583.0365789, 105.757208475, 0, 
                                                                                                                                                                                                                                                               2616.709167025, 2377.173263225, 75.694862075, 37.280938025, 34.3105211, 
                                                                                                                                                                                                                                                               0, 7804.60158, 4658.715853925, 266.33507635, 72.8675628, 0, 0, 
                                                                                                                                                                                                                                                               4168.6113629, 1692.9626727, 111.66667575, 0, 0, 0, 10032.36591165, 
                                                                                                                                                                                                                                                               2865.3011957, 1963.825942425, 69.45380285, 0, 0, 1785.1547225, 
                                                                                                                                                                                                                                                               999.668669225, 292.543231675, 254.033610875, 72.8069158, 0, 750.89466835, 
                                                                                                                                                                                                                                                               480.84163025, 331.199067375, 71.928040075, 37.61380175, 36.3270308, 
                                                                                                                                                                                                                                                               1083.2048175, 775.6721834, 681.0415219, 153.3758722, 0, 0, 854.93634015, 
                                                                                                                                                                                                                                                               825.524613275, 685.528790475, 141.46122475, 0, 0, 4818.203560075, 
                                                                                                                                                                                                                                                               833.174414525, 186.255344875, 35.771167625, 0, 0, 3932.2885342, 
                                                                                                                                                                                                                                                               2054.207895925, 524.696912425, 150.2404339, 0, 0, 2097.938841675, 
                                                                                                                                                                                                                                                               1579.990274075, 297.15061325, 35.444990125, 0, 0, 6577.38094975, 
                                                                                                                                                                                                                                                               5416.6666645, 505.95238075, 476.190476, 267.85714275, 178.5714285
                                                                                                                                                                                                                           )), class = "data.frame", row.names = c(NA, -108L), .Names = c("Crop", 
                                                                                                                                                                                                                                                                                          "Rot", "Rot.Herb", "species", "meandensity"))

# code 1
ggplot(weedseed3, aes(x=Rot.Herb, y=meandensity, fill=species))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set1", breaks=c("ABUTH","AMATA","CHEAL","SETSP","SOLPT","Others")) +
  theme_bw() +
  theme(panel.grid.major=element_blank()) +
  facet_grid(~Rot, scales = "free_x", space="free_x")+
  theme(legend.title=element_blank(), legend.text=element_text(size=20),legend.position="top")+
  xlab("\nTreatment") +
  theme(axis.title = element_text(size=24,face="bold", vjust=2), axis.text.x = element_text(size=20,angle = 90, hjust = 1, vjust=0.4)) +
  ylab("Seedbank density 2015 (seeds per meter sq)\n") +
  ylim(c(0,8000))+
  theme(axis.title = element_text(size=24,face="bold", vjust=2), axis.text.y = element_text(size=20, color="black"))

#  or following code2

ggplot(weedseed3, aes(x=Rot.Herb, y=meandensity, fill=species))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 6, name = "Set1"), 
                    breaks = c("ABUTH","AMATA","CHEAL","SETSP","SOLPT","Others")) +
  theme_bw() +
  theme(panel.grid.major=element_blank()) +
  facet_grid(~Rot, scales = "free_x", space="free_x")+
  theme(legend.title=element_blank(), legend.text=element_text(size=20),legend.position="top")+
  xlab("\nTreatment") +
  theme(axis.title = element_text(size=24,face="bold", vjust=2), axis.text.x = element_text(size=20,angle = 90, hjust = 1, vjust=0.4)) +
  ylab("Seedbank density 2015 (seeds per meter sq)\n") +
  ylim(c(0,8000))+
  theme(axis.title = element_text(size=24,face="bold", vjust=2), axis.text.y = element_text(size=20, color="black"))




