library(circlize)
library(biomaRt)
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(BioCircos)
library(grDevices)
library(ComplexHeatmap)

#Load your data
DatA <- readxl::read_excel('/media/charis/hdd/RNA-seq/DatA/res_DatA.xlsx')
#Download gene coordinates
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- DatA$GeneID
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description",
                              "chromosome_name", "start_position", "end_position"),
                values=genes,
                mart= mart
)
names(G_list) <- c("GeneID", "Entrez", "Symbol", "Description", "chr", "start","end") #
res_DatA <- merge(DatA, G_list, by = "GeneID")
res_DatA$chr_f = factor(res_DatA$chr,
                             levels=c('1','2','3','4','5','6',
                                      '7','8','9','10','11','12',
                                      '13','14','15','16','17','18',
                                      '19','20','21','22','X','Y','MT'))
res_DatA <- res_DatA[complete.cases(res_DatA$chr_f), ] #remove possible NAs

#create a bed-like dataframe for incorporation to the circlize package
bed_DatA = res_DatA[, c(15, 13, 14, 3 ,2)] #Hard coded (for reasons). Order should be: Chr, start, end and any value you want to plot against
bed_DatA$chr_f <- paste0("chr", bed_DatA$chr_f) #Needed for the circlize package.
bed_DatA$metap <- -log10(bed_DatA$metap) #Transform P to -log10(adj.p)

#Create two separate bed-like dataframes
bed_DatA_metafc <- bed_DatA[, c(1,2,3,4)]
bed_DatA_metafc$end <- bed_DatA_metafc$end+2000000
bed_DatA_metap <- bed_DatA[, c(1,2,3,5)]
bed_DatA_metap$metap[bed_DatA_metap$metap == "Inf"] <- 300 #Some P values may be pretty small and thus -log10(X) returns Inf.
bed_DatA_metap$end <- bed_DatA_metap$end+2000000 #Some genes are larger and thus ruin parts of the graph. Play around with this value to see what fits better.

#Same workflow for the second study
DatB <- readxl::read_excel('/media/charis/hdd/RNA-seq/DatB/res_DatB.xlsx')
res_DatB <- merge(DatB, G_list, by = "GeneID")
res_DatB$chr_f = factor(res_DatB$chr,
                             levels=c('1','2','3','4','5','6',
                                      '7','8','9','10','11','12',
                                      '13','14','15','16','17','18',
                                      '19','20','21','22','X','Y','MT'))
res_DatB <- res_DatB[complete.cases(res_DatB$chr_f), ]
bed_DatB = res_DatB[, c(15, 13, 14, 3 ,2)]
bed_DatB$chr_f <- paste0("chr", bed_DatB$chr_f)
bed_DatB$metap <- -log10(bed_DatB$metap)
bed_DatB_metafc <- bed_DatB[, c(1,2,3,4)]
bed_DatB_metafc$end <- bed_DatB_metafc$end+2000000
bed_DatB_metap <- bed_DatB[, c(1,2,3,5)]
bed_DatB_metap$metap[bed_DatB_metap$metap == "Inf"] <- 300
bed_DatB_metap$end <- bed_DatB_metap$end+2000000
bed_DatB_metap$chr_f <- paste0("DatB_", bed_DatB_metap$chr_f)
bed_DatA_metap$chr_f <- paste0("DatA_", bed_DatA_metap$chr_f)
bed_DatB_metafc$chr_f <- paste0("DatB_", bed_DatB_metafc$chr_f)
bed_DatA_metafc$chr_f <- paste0("DatA_", bed_DatA_metafc$chr_f)

#Bind both bed-like datasets from two situations 
df_metap <- rbind(bed_DatB_metap, bed_DatA_metap)
df_metafc <- rbind(bed_DatA_metafc, bed_DatB_metafc)

#create cytoband data and chromosome index
df_metafc <- rbind(bed_psoriasis_metafc, bed_Obesity_metafc)
chromosome.index = c(paste0("DatA_chr", c(1:22, "X", "Y")), 
                     rev(paste0("DatB_chr", c(1:22, "X", "Y")))) #Keep in mind that this refers to human genome.

#Graph initialization
circos.par(gap.after = c(rep(1, 23), 5, rep(1, 23), 5),
           start.degree = 90) #turn the graph around to look nicer.
circos.initializeWithIdeogram(df_metafc,#the above binding is a valid cytoband 
                              chromosome.index = chromosome.index,
                              plotType = NULL) #initialize cytoband
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
highlight.chromosome(paste0("DatA_chr", c(1:22, "X", "Y")), 
                     col = "red", track.index = 1) #Add Chr labels for DatA
highlight.chromosome(paste0("DatB_chr", c(1:22, "X", "Y")), 
                     col = "blue", track.index = 1) #Add Chr labels for DatB
circos.genomicTrack(df_metap, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value,
                       col = ifelse(value[[1]] > 10, "#FF5733", "black"),
                       cex = 0.5, ...)
                    }) #Add P value scatter plot
circos.genomicTrack(df_metafc,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "red", "cyan"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 5, col = "white")
                    }) #Add logFC barplot

#Add some labels
text(-0.9, -0.8, "DatB\nRNA-seq results") #Description of the left circoplot
text(0.9, 0.8, "DatA\nRNA-seq results") #Description of the right circoplot
text(0, 0, expression(paste("Gene/Chromosome adjusted ",italic("P")* " and log"[2]*"FC results"))) #Center your title (or add any image you might want)
