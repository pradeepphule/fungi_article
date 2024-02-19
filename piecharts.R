# Simple Pie Chart
slices <- c(10, 12,4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie(slices, labels = lbls, main="Pie Chart of Countries")

# Pie Chart with Percentages
slices <- c(10, 12, 4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries") 

# 3D Exploded Pie Chart
library(plotrix)
slices <- c(10, 12, 4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie3D(slices,labels=lbls,explode=0.2,
      main="Pie Chart of Countries ") #value of explode can be increase separation of pie
df <- data.frame(pcategory = c('Production','Facilities','Labor','Licenses','Taxes',
                               'Legal','Insurance'), pcost = c(24.8,20.98,17.48,
                                                             12.59,10.49,8.39,5.59))
library(ggvis)
library(dplyr)
library(scales)
df <- data.frame(pcategory = c('Production','Facilities','Labor','Licenses','Taxes',
                               'Legal','Insurance'), pcost = c(24.8,20.98,17.48,
                                                               12.59,10.49,8.39,5.59))
#Sort from highest to lowest, levels are also reordered
df <- transform(df, pcategory=reorder(pcategory, -pcost) ) 
#Draw bar chart
df %>% ggvis(~pcategory,~pcost, fill = "blue") %>% layer_bars() %>%
  add_axis("x", title = "Project Category") %>%
  add_axis("y", title = "Project Cost") %>%
  add_axis("x", orient = "top", ticks = 0, title = "Project Cost Breakdown", #title
           properties = axis_props(
             axis = list(stroke = "white"),
             labels = list(fontSize = 0))) %>%
  hide_legend("fill") #Hide legends


pieval<-c(2,4,6,8)
pielabels<-
  c("We hate\n pies","We oppose\n  pies","We don't\n  care","We just love pies")
# grab the radial positions of the labels
lp<-pie3D(pieval,radius=0.9,labels=pielabels,explode=0.1,main="3D PIE OPINIONS")
# lengthen the last label and move it to the left
pielabels[4]<-"We cannot survive without our pies"
lp[4]<-4.8
# specify some new colors
pie3D(pieval,radius=0.9,labels=pielabels,explode=0.1,main="3D PIE OPINIONS",
      col=c("brown","#ddaa00","pink","#dd00dd"),labelpos=lp)



##############################################################

#below is rough code for exerciese for idea only
-----------------------------------------------------------------------------#
  #-----------------------------------------------------------------------------#
  #                 R in Action - Basic Graphs 
  #             - Bar Plots
  #             - Pie charts
  #             - Fan plots
  #             - Histogram
  #             - Kernel Density Plot
  #             - Box plot
  #             - Dot plot
  #-----------------------------------------------------------------------------#
  #-----------------------------------------------------------------------------#


install.packages(c('vcd','plotrix','sm','vioplot'))

# pause after each graph
par(ask = TRUE)

# save original graphic settings
opar <- par(no.readonly = TRUE)

# Load vcd package in order to get teh Arthritis dataset
library(vcd)


#############################################
###    ----------- Bar Plots ----------   ###
#############################################

####################
# 1. simple bar plot
####################

# Simple bar plot
counts <- table(Arthritis$Improved) # Get cell counts first
counts

barplot(counts, main = "Simple Bar Plot", xlab = "Improvement", 
        ylab = "Frequency") # Plot

# Horizontal bar plot




##################################
# 2. Stacked and groupde bar plots
################################## 

counts <- table(Arthritis$Improved, Arthritis$Treatment) # Get cell counts first
counts

# stacked barplot
barplot(counts, main = "Stacked Bar Plot", xlab = "Treatment", 
        ylab = "Frequency", col = c("red", "yellow", "green"), 
        legend = rownames(counts))

# grouped barplot
barplot(counts, main = "Grouped Bar Plot", xlab = "Treatment", 
        ylab = "Frequency", col = c("red", "yellow", "green"), 
        legend = rownames(counts), 
        beside = TRUE)


##################################
# 3. Mean bar plots
################################## 

states <- data.frame(state.region, state.x77)
means <- aggregate(states$Illiteracy, 
                   by = list(state.region), 
                   FUN = mean)
means

means <- means[order(means$x), ]
means

barplot(means$x, names.arg = means$Group.1)
title("Mean Illiteracy Rate")

##################################
# 4. Fitting labels in bar plots
##################################


par(mar = c(5, 8, 4, 2))
par(las = 2) # labels are parallel (=0) or perpendicular(=2) to axis
counts <- table(Arthritis$Improved)

barplot(counts, main = "Treatment Outcome", horiz = TRUE, 
        cex.names = 0.8, names.arg = c("No Improvement", 
                                       "Some Improvement", "Marked Improvement"))

##################################
# 4. Spinograms use for publication
##################################
# Spinograms is a stacked bar plot is rescaled so that the height of
#each bar is 1 and the segment heights represent proportions

library(vcd)
attach(Arthritis)
counts <- table(Treatment, Improved)
spine(counts, main = "Spinogram Example")
detach(Arthritis)


#############################################
###    ----------- Pie Charts ---------   ###
#############################################

# Pie chart with Frequency
par(mfrow = c(2, 2))
slices <- c(10, 12, 4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")

pie(slices, labels = lbls, main = "Simple Pie Chart")

# Percentage as label
pct <- round(slices/sum(slices) * 100)
lbls2 <- paste(lbls, " ", pct, "%", sep = "")
pie(slices, labels = lbls2, col = rainbow(length(lbls)), 
    main = "Pie Chart with Percentages")

# 3D Pie
library(plotrix)
pie3D(slices, labels = lbls, explode = 0.5, main = "3D Pie Chart ")

# Pie chart from a table
mytable <- table(state.region)
mytable
lbls <- paste(names(mytable), "\n", mytable, sep = "")
pie(mytable, labels = lbls, 
    main = "Pie Chart from a Table\n (with sample sizes)")


par(opar)


#############################################
###    ----------- Fan Plots  ---------   ### use this for publication
#############################################

# Pie charts make it difficult to compare the values of the slices (unless the values are appended to the labels).

library(plotrix)
slices <- c(10, 12, 4, 16, 8)
lbls <- c("US", "UK", "Australia", "Germany", "France")
fan.plot(slices, labels = lbls, col = rainbow(length(lbls)) , main = "Fan Plot")



#############################################
###    ----------- Histogram  ---------   ###
#############################################

# Histograms

par(mfrow = c(2, 2))

hist(mtcars$mpg)

# Specify # of bins
hist(mtcars$mpg, breaks = 12, col = "red", 
     xlab = "Miles Per Gallon", 
     main = "Colored histogram with 12 bins")

# Histogram with density curve and rug plot overlay
hist(mtcars$mpg, freq = FALSE, breaks = 12, col = "red", 
     xlab = "Miles Per Gallon", 
     main = "Histogram, rug plot, density curve")
rug(jitter(mtcars$mpg))
lines(density(mtcars$mpg), col = "blue", lwd = 2)
a <- density(mtcars$mpg)

# Histogram with Superimposed Normal Curve 
x <- mtcars$mpg
h <- hist(x, breaks = 12, col = "red", 
          xlab = "Miles Per Gallon", 
          main = "Histogram with normal curve and box")
xfit <- seq(min(x), max(x), length = 40)
yfit <- dnorm(xfit, mean = mean(x), sd = sd(x))
yfit <- yfit * diff(h$mids[1:2]) * length(x)
lines(xfit, yfit, col = "blue", lwd = 2)
box()

# restore original graphic parameters
par(opar)

#############################################
###   ------ Kernel density plot  -----   ###
#############################################

# 1. Kernel density plots
par(mfrow = c(2, 1))
d <- density(mtcars$mpg)

plot(d)

d <- density(mtcars$mpg)
plot(d, main = "Kernel Density of Miles Per Gallon")
polygon(d, col = "red", border = "blue")
rug(mtcars$mpg, col = "brown")

par(opar)

# 2. Comparing multiple kernel density plots

par(lwd = 2)
library(sm)
attach(mtcars)

mtcars
cyl.f <- factor(cyl, levels = c(4, 6, 8), 
                labels = c("4 cylinder", "6 cylinder", "8 cylinder"))

sm.density.compare(mpg, cyl, xlab = "Miles Per Gallon")
title(main = "MPG Distribution by Car Cylinders")

colfill <- c(2:(1 + length(levels(cyl.f))))
cat("Use mouse to place legend...", "\n\n")
legend(locator(1), levels(cyl.f), fill = colfill)
detach(mtcars)
par(lwd = 1)


#############################################
###   ------------ Box Plot  ----------   ###
#############################################

####################
# 1. Simple box plot
####################

boxplot(mtcars$mpg, main="Box plot", ylab="Miles per Gallon")

boxplot.stats(mtcars$mpg)


##########################
# 2. Using parallel box plot to compare groups
############################
mtcars
boxplot(mpg ~ cyl, data = mtcars, 
        main = "Car Milage Data", 
        xlab = "Number of Cylinders", 
        ylab = "Miles Per Gallon")

# Adding the option varwidth=TRUE will make the box plot widths proportional
#to the square root of their sample sizes.
boxplot(mpg ~ cyl, data = mtcars, notch = TRUE, 
        varwidth = TRUE, col = "red", 
        main = "Car Mileage Data", 
        xlab = "Number of Cylinders", 
        ylab = "Miles Per Gallon")  


##############################
# 3. Box plots for two crossed factors
################################

mtcars$cyl.f <- factor(mtcars$cyl, levels = c(4, 6, 8), labels = c("4", "6", "8"))

mtcars$am.f <- factor(mtcars$am, levels = c(0, 1), 
                      labels = c("auto", "standard"))

boxplot(mpg ~ am.f* cyl.f, data = mtcars, 
        varwidth = TRUE, col = c("gold", "darkgreen"), 
        main = "MPG Distribution by Auto Type", 
        xlab = "Auto Type")



#############################################
###   ------------ Dot Plot  ----------   ###
#############################################

# 1. dot plot
dotchart(mtcars$mpg, labels = row.names(mtcars), 
         cex = 0.7, 
         main = "Gas Milage for Car Models", 
         xlab = "Miles Per Gallon")

# 2. sorted colored grouped dot chart

x <- mtcars[order(mtcars$mpg), ]
x$cyl <- factor(x$cyl)
x$color[x$cyl == 4] <- "red"
x$color[x$cyl == 6] <- "blue"
x$color[x$cyl == 8] <- "darkgreen"

dotchart(x$mpg, labels = row.names(x), cex = 0.7, 
         pch = 19, groups = x$cyl, 
         gcolor = "black", color = x$color,
         main = "Gas Milage for Car Models\ngrouped by cylinder", 
         xlab = "Miles Per Gallon")

