library(ggplot2)
library(plyr)
library(randomcoloR)
library(readr)
rm(list = ls(all = TRUE))
data <- read_csv("~/Desktop/biological_plots/interproscan_repeats_top20", 
                 col_names = TRUE)

data
View(data)
obj <- ggplot(data,aes(x = Fusarium, y = Seqs, fill = `IPS Repeats`,label = Seqs)) +
  geom_bar(colour="white", stat="identity" )+
  geom_text(size = 3, position = position_stack(vjust = 0.5))
n <- 20
palette <- distinctColorPalette(n)

obj + scale_fill_manual(values = palette, 
                        breaks = data$`IPS Repeats`)
