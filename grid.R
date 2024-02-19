set.seed(123)
mat = matrix(sample(16, 16), 4, byrow = TRUE)
rownames(mat) = paste0("R", 1:4)
colnames(mat) = paste0("C", 1:4)
mat
df = reshape2::melt(mat)
df
grid.col = c("#00000040", "red", "#00FF0040", "#0000FF40", "orange", "pink", "yellow", "grey")
grid.col
names(grid.col) = c(rownames(mat), colnames(mat))
grid.col
library(circlize)
chordDiagram(df, grid.col = grid.col)