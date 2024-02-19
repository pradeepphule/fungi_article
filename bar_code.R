library(gplots)
cil <- data_final3 * 0.85  
ciu <- data_final3 * 1.15  
barplot2(data_final3,                             # Data (bar heights) to plot  
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
legend("topleft",                                # Add a legend to the plot  
legend=c("With salt medium","Without salt medium","Marine bacterium culture","Mycobacterium vaccae","E.coli","Bacillus subtillus"),
 fill=c("lightblue", "mistyrose",
 "lightcyan", "lavender","cornsilk","grey"), bty = "n")                  # Fill for boxes of the legend  
