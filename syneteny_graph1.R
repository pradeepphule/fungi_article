library(genoPlotR)
df1 <- data.frame(name=c
                  ("ctg55_orf40","ctg55_orf41","ctg55_orf42","ctg55_orf43","ctg55_orf44","ctg55_orf45","ctg55_orf46","ctg55_orf47","ctg55_orf48","ctg55_orf49","ctg55_orf50","ctg55_orf51","ctg55_orf52","ctg55_orf53")
                  ,start=c(2258,3644,4850,7438,16284,18251,20001,27591,29435,34043,42018,43066,44289,47636)
                  ,end=c(3478,4450,6920,11316,17721,19623,26872,28511,33160,34970,42901,43967,46837,48308)
                  ,strand=c(1,-1,1,-1,1,-1,1,1,-1,1,1,-1,1,-1)
                  ,col=c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"))

dna_seg1 <- dna_seg(df1)

df2 <- data.frame(name=c
                  ("CZR48045","CZR48046","CZR48047","CZR48048","CZR48049","CZR48050","CZR48051","CZR48052","CZR48053","CZR48054","CZR48055","CZR48056","CZR48057","CZR48058","CZR48059","CZR48060","CZR48061","CZR48062")
                  ,start=c(133811,135606,137732,140581,142593,144422,149121,151173,153396,154244,156870,159465,162679,164011,166660,168825,172836,175526)
                  ,end=c(135385,136719,138352,141887,144050,146747,151072,153002,154211,156516,158211,162609,163451,166069,168315,171593,174973,176615)
                  ,strand=c(1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,-1)
                  ,col=c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"))

dna_seg2 <- dna_seg(df2)

df3 <- data.frame(name=c
                  ("CCT73160","CCT73161","CCT73162","CCT73163","CCT73164","CCT73165","CCT73166","CCT73167","CCT73168","CCT73169","CCT73170","CCT73171","CCT73172","CCT73173","CCT73174","CCT73175")
                  ,start=c(99483,100446,103979,114408,116825,118877,121103,121952,124583,127362,130371,131533,134181,136346,140358,143049)
                  ,end=c(100103,102548,105285,115865,118776,120708,121918,124031,125920,130301,131140,133522,135834,139116,142492,144138)
                  ,strand=c(-1,-1,1,1,1,-1,1,-1,-1,1,-1,-1,-1,-1,1,-1)
                  ,col=c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"))

dna_seg3 <- dna_seg(df3)

df4 <- data.frame(name=c
                  ("XP_018258548","XP_018258549","XP_018258550","XP_018258551","XP_018258552","XP_018258553","XP_018258554","XP_018258555","XP_018258556")
                  ,start=c(86389,88994,93747,95251,98460,100649,101855,103827,106295)
                  ,end=c(87772,93150,95117,97225,100630,101462,103684,105789,106390)
                  ,strand=c(1,1,1,-1,1,-1,1,-1,1)
                  ,col=c("blue","blue","blue","blue","blue","blue","blue","blue","blue"))

dna_seg4 <- dna_seg(df4)

dna_segs <- list(dna_seg1, dna_seg2, dna_seg3, dna_seg4)
names <- c("Huey", "Dewey", "Louie","Pradeep")
names(dna_segs) <- names
#("ctg55_orf40","ctg55_orf41","ctg55_orf42")
cmp1 <- data.frame(start1=c(2258,3644,4850),
                   end1=c(3478,4450,6920),
                   #("CZR48054","CZR48053","CZR48052")
                   start2=c(154244,153396,151173),
                   end2=c(156516,154211,153002),
                   col=c("#67000D", "#08306B", "#08306B"))

comparison1 <- comparison(cmp1)

#("CZR48054","CZR48053","CZR48052")
cmp2 <- data.frame(start1=c(154244,153396,151173),
                   end1=c(156516,154211,153002),
                   #("CCT73167","CCT73166","CCT73165")
                   start2=c(121952,121103,118877),
                   end2=c(124031,121918,120708),
                   col=c("#67000D", "#08306B", "#08306B"))

comparison2 <- comparison(cmp2)


#("CCT73166","CCT73165")
cmp3 <- data.frame(start1=c(121103,118877),
                    end1=c(121918,120708),
                    #("XP_018258553","XP_018258554)
                    start2=c(100649,101855),
                    end2=c(101462,103684),
                    col=c("#67000D", "#08306B"))
                   
                   comparison3 <- comparison(cmp3)
                   
                   
                   comparisons <- list(comparison1, comparison2, comparison3)
                   
                   plot_gene_map(dna_segs=dna_segs, comparisons=comparisons, main="Comparison of Huey, Dewey and Louie")
                   