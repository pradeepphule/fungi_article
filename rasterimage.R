set.seed(123)
mat = matrix(rnorm(10000*100), ncol = 100)
library(ComplexHeatmap)
pdf("heatmap.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()
pdf("heatmap_raster_by_png.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE, raster_device = "png")
dev.off()

pdf("heatmap_raster_by_jpeg.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE, raster_device = "jpeg")
dev.off()

pdf("heatmap_raster_by_tiff.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE, raster_device = "tiff")
dev.off()

pdf("heatmap_raster_by_CairoPNG.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE, raster_device = "CairoPNG")
dev.off()

pdf("heatmap_raster_by_CairoJPEG.pdf", width = 8, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE, raster_device = "CairoJPEG")
dev.off()

all_files = c("heatmap.pdf", "heatmap_raster_by_png.pdf", 
              "heatmap_raster_by_jpeg.pdf", "heatmap_raster_by_tiff.pdf",
              "heatmap_raster_by_CairoPNG.pdf", "heatmap_raster_by_CairoJPEG.pdf")
fs = file.size(all_files)
names(fs) = all_files
sapply(fs, function(x) paste(round(x/1024), "KB"))
