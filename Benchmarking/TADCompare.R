library(TADCompare)
library(dplyr)
library(writexl)

resoultion = 40000
# read both file from HiCPro
mat_pts <- read.table('PTs_40000_iced.matrix')
bed_pts <- read.table('PTs_40000_abs.bed')

mat_rts <- read.table('RTs_40000_iced.matrix')
bed_rts <- read.table('RTs_40000_abs.bed')

# Convert to modified bed format
sparse_mats_pts <- HiCcompare::hicpro2bedpe(mat_pts,bed_pts)
sparse_mats_rts <- HiCcompare::hicpro2bedpe(mat_rts,bed_rts)

saveRDS(sparse_mats_pts, 'sparse_mats_pts.rds')
saveRDS(sparse_mats_rts, 'sparse_mats_rts.rds')
# Remove empty matrices if necessary
# sparse_mats$cis = sparse_mats$cis[sapply(sparse_mats, nrow) != 0]

sparse_mats_pts <- readRDS('sparse_mats_pts.rds')
sparse_mats_rts <- readRDS('sparse_mats_rts.rds')
# load TADs
pts_tads <- read.table('/data/kfang/BRCA_HiC_Tissue/revision/GISTA_benchmark/TopDom/PTs/PTs_40000_win5_total.txt', header = T)
rts_tads <- read.table('/data/kfang/BRCA_HiC_Tissue/revision/GISTA_benchmark/TopDom/RTs/RTs_40000_win5_total.txt', header = T)

# only keep tads
pts_tads <- pts_tads[pts_tads$tag == 'domain', ]
rts_tads <- rts_tads[rts_tads$tag == 'domain', ]
pts_tads <- pts_tads[,c('chr','from.coord','to.coord')]
rts_tads <- pts_tads[,c('chr','from.coord','to.coord')]
colnames(pts_tads) <- c('chr','start','end')
colnames(rts_tads) <- c('chr','start','end')

# Go through all pairwise chromosomes and run TADCompare
sparse_tads = lapply(1:length(sparse_mats_pts$cis), function(z) {
  x <- sparse_mats_pts$cis[[z]]
  y <- sparse_mats_rts$cis[[z]]
  
  #Pull out chromosome
  chr <- x[, 1][1]
  #Subset to make three column matrix
  x <- x[, c(2, 5, 7)]
  y <- y[, c(2, 5, 7)]
  #Run SpectralTAD
  combined_bed <- list(pts_tads[pts_tads$chr == chr,], rts_tads[rts_tads$chr == chr,])
  comp <- TADCompare(x, y, resolution = 40000, pre_tads = combined_bed)
  return(list(comp, chr))
})

# Pull out differential TAD results
diff_res <- lapply(sparse_tads, function(x) x[[1]])
# Pull out chromosomes
chr <- lapply(sparse_tads, function(x) x[[2]])
# Name list by corresponding chr
names(diff_res) <- chr
chromosomes <- names(diff_res)
# Combine TAD_Frame data for all chromosomes in diff_res, handling empty data frames
combined_TAD_Frame <- do.call(rbind, lapply(chromosomes, function(chr) {
  # Access TAD_Frame for the current chromosome
  tad_frame <- diff_res[[chr]][['TAD_Frame']]
  # Check if the TAD_Frame is not empty
  if (nrow(tad_frame) > 0) {
    # Add a column for chromosome
    tad_frame$chr <- chr
    # Reorder columns to move 'chr' to the first position
    tad_frame <- tad_frame[, c("chr", setdiff(names(tad_frame), "chr"))]
  }
  # Return the TAD_Frame (empty or populated)
  return(tad_frame)
}))

combined_Boundary_Score <- do.call(rbind, lapply(chromosomes, function(chr) {
  # Access TAD_Frame for the current chromosome
  tad_frame <- diff_res[[chr]][['Boundary_Scores']]
  # Check if the TAD_Frame is not empty
  if (nrow(tad_frame) > 0) {
    # Add a column for chromosome
    tad_frame$chr <- chr
    # Reorder columns to move 'chr' to the first position
    tad_frame <- tad_frame[, c("chr", setdiff(names(tad_frame), "chr"))]
  }
  # Return the TAD_Frame (empty or populated)
  return(tad_frame)
}))

# save to excel
sheets <- list('TAD_Frame'=combined_TAD_Frame, 'Boundary_Score'=combined_Boundary_Score)
write_xlsx(sheets,'PTsRTs_TADCompare_results.xlsx')