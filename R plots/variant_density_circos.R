library(circlize)
library(tidyverse)


par(mar = c(1, 1, 3, 8))  # Bottom, Left, Top, Right margins

circos.clear()

# Read data
genome <- read.table("genome_structure.txt", header = FALSE)
colnames(genome) <- c("chr", "start", "end")

variants <- read.table("variants_hq.bed", header = FALSE)
colnames(variants) <- c("chr", "start", "end", "type", "qual")

coverage <- read.table("coverage_windows.txt", header = FALSE)
colnames(coverage) <- c("chr", "start", "end", "depth")

variant_density <- read.table("variant_density.txt", header = FALSE)
colnames(variant_density) <- c("chr", "start", "end", "count")

print(paste("Variant density windows:", nrow(variant_density)))
print(paste("Max variants per window:", max(variant_density$count)))

# Robust scaling
MAX_SPIKE <- quantile(variant_density$count, 0.99)
MAX_COV <- 200

# Setup
circos.par(
  start.degree = 90,
  gap.degree = 10,
  cell.padding = c(0, 0, 0, 0)
)



circos.genomicInitialize(genome, plotType = NULL)

# Track 1: Labels
circos.track(ylim = c(0, 1), 
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               circos.text(CELL_META$xcenter, 0.5, chr, 
                           cex = 1, facing = "bending.inside",
                           niceFacing = TRUE)
             }, 
             track.height = 0.1, bg.border = NA)

# Track 2: Coverage
circos.genomicTrack(coverage,
                    panel.fun = function(region, value, ...) {
                      colors <- colorRampPalette(c("white", "yellow", "red"))(100)
                      col_idx <- pmax(1, pmin(100, round((value[[1]]/MAX_COV) * 100)))
                      circos.genomicRect(region, value, col = colors[col_idx], border = NA)
                    },
                    track.height = 0.15)

# Read SNP and indel densities
snp_density <- read.table("snp_density.txt", header = FALSE)
colnames(snp_density) <- c("chr", "start", "end", "count")

indel_density <- read.table("indel_density.txt", header = FALSE)
colnames(indel_density) <- c("chr", "start", "end", "count")

# Robust scaling (use max from both)
MAX_SPIKE <- max(quantile(snp_density$count, 0.99), 
                 quantile(indel_density$count, 0.99))

# Track 3: Overlaid SNP (blue) and Indel (red) spikes
circos.trackPlotRegion(
  ylim = c(0, MAX_SPIKE),
  track.height = 0.25,
  bg.col = "#F7F7F7",
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    
    # Draw SNPs (blue) first
    snp_data <- snp_density[snp_density$chr == chr, ]
    if(nrow(snp_data) > 0) {
      xmid <- (snp_data$start + snp_data$end) / 2
      heights <- pmin(snp_data$count, MAX_SPIKE)
      circos.segments(x0 = xmid, y0 = 0, x1 = xmid, y1 = heights,
                      col = "#0066FF", lwd = 1)
      circos.points(x = xmid, y = heights, pch = 16, cex = 0.3, col = "#0066FF")
    }
    
    # Draw indels (red) on top
    indel_data <- indel_density[indel_density$chr == chr, ]
    if(nrow(indel_data) > 0) {
      xmid <- (indel_data$start + indel_data$end) / 2
      heights <- pmin(indel_data$count, MAX_SPIKE)
      circos.segments(x0 = xmid, y0 = 0, x1 = xmid, y1 = heights,
                      col = "#FF3333", lwd = 1.2)
      circos.points(x = xmid, y = heights, pch = 17, cex = 0.4, col = "#FF3333")
    }
  }
)

title("Salmonella enterica Genome - Variant Density (Q>=20)")


legend("right", 
       legend = c("SNPs (5kb bins)", "Indels (5kb bins)", "Low coverage", "High coverage"),
       col = c("#0066FF", "#FF3333", "yellow", "red"),
       lwd = c(2, 2, NA, NA),
       pch = c(16, 17, 15, 15),
       pt.cex = c(1, 1, 2, 2),
       cex = 0.9,
       bty = "n")

