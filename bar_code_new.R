library(gplots)
#Clean up workspace - i.e. delete variable created by the graphics demo
rm(list = ls(all = TRUE))
working_dir = "/home/pradeep/fusarium/analysis/0Growth_data" 
setwd(working_dir)
data_final = read.csv("data_final.csv")
upper_bound = read.csv("upper_bound.csv")
lower_bound = read.csv("lower_bound.csv")
data_final3 = as.matrix(data_final1)
cil <- t(as.matrix(upper_bound))  
ciu <- t(as.matrix(lower_bound))  
barplot2(t(as.matrix(data_final)),                           # Data (bar heights) to plot  
    beside=TRUE,                            # Plot the bars beside one another; default is to plot stacked bars  
    space=c(0.2,0.8),                       # Amount of space between i) bars within a group, ii) bars between   
    names.arg=c("1", "2", "3", "4", "5" ,"6" ,"7"),    #Names for the bars  
    col=c("lightblue", "mistyrose" ,"lightcyan","lavender","cornsilk","grey"),                   # Color of the bars  
    border="black",                         # Color of the bar borders  
    main=c("Growth of Fusarium pp1-1"),      # Main title for the plot  
    xlab="time(days)",                       # X-axis label  
    ylab="growth(mm)", ylim = c(0,50),las = 1,     # cex.axis = 1, cex.names = 1,             # Y-axis label  
    font.lab=2,                            # Font to use for the axis labels: 1=plain text, 2=bold,
    plot.ci=TRUE,                           # Plot confidence intervals  
    ci.l=cil,                               # Lower values for the confidence interval  
    ci.u=ciu,                               # Upper values for the confidence interval  
    plot.grid=FALSE)                         # Plot a grid  
axis(side = 2, at = seq(5,45,5) , labels = FALSE)
legend("topleft",                                
legend=c("With salt medium","Without salt medium","Marine bacterium culture","Mycobacterium vaccae","E.coli","Bacillus subtillus"),
 fill=c("lightblue", "mistyrose","lightcyan", "lavender","cornsilk","grey"), bty = "n")     # Add a legend to the plot               # Fill for boxes of the legend  
