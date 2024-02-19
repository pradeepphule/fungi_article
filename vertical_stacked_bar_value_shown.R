library(ggplot2)
library(plyr)
library(randomcoloR)
library(readr)
rm(list = ls(all = TRUE))
data <- read_csv("~/Desktop/biological_plots/interproscan_family_top20", 
                 col_names = TRUE)

data
View(data)
obj <- ggplot(data,aes(x = Fusarium, y = Seqs, fill = `IPS Family`,label = Seqs)) +
  geom_bar(colour="white", stat="identity" )+
  geom_text(size = 3, position = position_stack(vjust = 0.5))
n <- 20
palette <- distinctColorPalette(n)

obj + scale_fill_manual(values = palette, 
                        breaks = data$`IPS Family`)







 
 
 
 
 
 
#Year      <- c(rep(c("2006-07", "2007-08", "2008-09", "2009-10"), each = 4))
#Category  <- c(rep(c("A", "B", "C", "D"), times = 4))
#Frequency <- c(168, 259, 226, 340, 216, 431, 319, 368, 423, 645, 234, 685, 166, 467, 274, 251)
#Data      <- data.frame(Year, Category, Frequency)
#Data
#data<-t(data)
#colnames(data) <- data[1, ]

#data<-as.data.frame(data)
#data
#ggplot(data, aes(x = Fusarium, y = Seqs, fill = `IPS Family`, label = Seqs)) +
#ggplot(data, aes(x = Year, y = Frequency, fill = Category, label = Frequency)) + # to make x and y axis
#  geom_bar(colour="white",stat = "identity") +                              # to show the bars
#geom_text(size = 3, position = position_stack(vjust = 0.5))# to show the values 

# allcolor<-c("A" = "black", "B" = "orange", "C" = "blue",
#             "D" = "yellow", "E" = "thistle4", "F" = "thistle1",
#             "G" = "yellow3", "H" = "tomato1", "I" = "tan3",
 #            "J" = "yellowgreen", "K" = "tomato4", "L" = "tan",
#             "M" = "wheat4", "N" = "turquoise4", "O" = "steelblue",
 #            "P" = "wheat3", "Q" = "violet", "R" = "springgreen4",
#             "S" = "wheat1", "T" = "violetred3")
# 
# obj + scale_fill_manual(values = allcolor,breaks = data$`IPS Family`[1])
# obj + scale_fill_manual(values = RColorBrewer::brewer.pal(n = 9, name = "Set1"), 
#                         breaks = data$`IPS Family`)
 
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# obj + scale_fill_manual(values = color, 
 #                        breaks = data$`IPS Family`)