library(genoPlotR)

#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))
###########################################################
data(three_genes)
comparisons[[1]]$col <- apply_color_scheme(c(0.6, 0.4, 0.5), "grey")
names <- c("Huey", "Dewey", "Louie")
names(dna_segs) <- names
tree <- newick2phylog("(((Huey:4.2,Dewey:3.9):3.1,Louie:7.3):1);")
mid_pos <- middle(dna_segs[[1]])
xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(1850, 2800))
annot <- annotation(x1=c(mid_pos[1], dna_segs[[1]]$end[2]),
                    x2=c(NA, dna_segs[[1]]$end[3]),
                    text=c(dna_segs[[1]]$name[1], "region1"),
                    rot=c(30, 0), col=c("blue", "black"))
plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
                annotations=annot, annotation_height=1.3,
                tree=tree, tree_width=2,
                xlims=xlims,
                main="Comparison of Huey, Dewey and Louie")

#######################################################################
df1 <- data.frame(name=c("feat1", "feat2", "feat3"),
                  start=c(2, 1000, 1050),
                  end=c(600, 800, 1345),
                  strand=c(-1, -1, 1),
                  col=c("blue", "black", "red"))
dna_seg1 <- dna_seg(df1)
df2 <- data.frame(name=c("feat4", "feat5", "feat6"),
                  start=c(50, 800, 1200),
                  end=c(900, 1100, 1322),
                  strand=c(-1, 1, 1),
                  col=c("blue", "black", "red"))
dna_seg2 <- dna_seg(df2)
df3 <- data.frame(name=c("feat7", "feat8", "feat9"),
                  start=c(1899, 2108, 2803),
                  end=c(2034, 2732, 3620),
                  strand=c(-1, -1, 1),
                  col=rep("blue", 3))
dna_seg3 <- dna_seg(df3)
dna_segs <- list(dna_seg1, dna_seg2, dna_seg3)

df4 <- data.frame(start1=dna_seg1$start,
                    end1=dna_seg1$end,
                    start2=dna_seg2$start,
                    end2=dna_seg2$end)
 comparison1 <- comparison(df4)
 df5 <- data.frame(start1=c(50, 800),
                      end1=c(500, 1100),
                      start2=c(1899, 2732),
                      end2=c(2034, 2508),
                      col=c("#67000D", "#08306B"))
 comparison2 <- comparison(df5)
 comparisons <- list(comparison1, comparison2)
 
 plot_gene_map(dna_segs=dna_segs, comparisons=comparisons)
 
 
