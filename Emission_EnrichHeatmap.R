library(genomation)
library(pheatmap)
library('biomaRt')
ensembl <- useMart('ensembl')  # using ensembl database data
mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
library(ensembldb)
library(EnsDb.Mmusculus.v79)
names(listTables(EnsDb.Hsapiens.v75))
edb <- EnsDb.Mmusculus.v79  # abbreviate
ensembl_genes <- genes(edb)

### =====emission grid heatmap================================================================
library(ggplot2)
emission_df <- read.delim('./emissions_6.txt')
emission_mat <- as.matrix(emission_df[1:6, 2:6])
### reorder emission_mat columns, then transpose matrix
emission_mat <- t(emission_mat[ , c(4,5,3,1,2)])
### reorder rows 
emission_mat <- emission_mat[c(4,5,3,1,2), ]
library("RColorBrewer")
library(EnrichedHeatmap)
col_fun <- colorRampPalette(brewer.pal(10, "GnBu"))(256)

Heatmap(emission_mat, name = "mat", col = col_fun, cluster_rows=F, cluster_columns=F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.4f", emission_mat[i, j]), x, y, 
                    gp = gpar(fontsize = 12))})
### =========================================================================================
### seperate in-house 6 states
# add strand as '*' into column 4 before readBed
ES_file <- readBed('./ChromHMM_inhouse/Mouse_liver_6_segments.bed')
head(ES_file)
ES1 <- ES_file[grepl('E1', ES_file$score), ]
ES2 <- ES_file[ES_file$score == 'E2', ]
ES3 <- ES_file[ES_file$score == 'E3', ]
ES4 <- ES_file[ES_file$score == 'E4', ]
ES5 <- ES_file[ES_file$score == 'E5', ]
ES6 <- ES_file[ES_file$score == 'E6', ]
### put ES together
ES <- c(ES1, ES2, ES3, ES4, ES5, ES6)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library("org.Mm.eg.db")
library(ChIPseeker)

Ann <- as.data.frame(annotatePeak(ES, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                  annoDb = "org.Mm.eg.db"))
head(Ann)
library(dplyr)
### remove duplicates, keep the state closest to tss
## GL: gene list from state annotation
GL <- as.data.frame(Ann %>% group_by(ENSEMBL) %>% dplyr::filter(distanceToTSS == min(abs(distanceToTSS))))
### multiple transcripts' distanceToTSS could be all 0, then take the most upstream one
GL <- as.data.frame(GL %>% group_by(ENSEMBL) %>% dplyr::filter(geneStart == min(geneStart)))
GL[duplicated(GL$ENSEMBL), ]
### there are still some duplicates, which must be in the same state
GL <- GL[!duplicated(GL$ENSEMBL), ]
### add expression to GL
### expression
wt_fpkm <- read.csv('../Mouse_RNAseq/FPKMs/Mfuzz/WT_FPKM_mean.csv', sep='\t', header=T)
fpkm_bl <- data.frame(ensembl=rownames(wt_fpkm), bl_exp=wt_fpkm[,1])

### Sig_DEG_afterPH
library(xlsx)
sig_genes <- c()
for ( i in 1:8){ 
  dat <- read.xlsx('../Regeneration/Sig_DEG_afterPH_eachT.xlsx', sheetIndex = i)
  sig_genes <- c(sig_genes, as.character(dat$ensembl))
  sig_genes <- unique(sig_genes) # 2179 sig genes in 7 comparision
}

### subset ensemble_genes with ensembl id in expression spreadsheet
GL_exp <- merge(GL, fpkm_bl, by.x='ENSEMBL', by.y='ensembl', all.x=T)
### remove those are not detected 
GL_exp <- GL_exp[!is.na(GL_exp$bl_exp), ]
### make GRagnes of GL (gene list) for plotting Chip-seq signals
sub_genes <-  subset(ensembl_genes, gene_id %in% GL_exp$ENSEMBL)
gr_df <- as.data.frame(sub_genes)
## 1: ENSEMBL, 7: score(state), 20:bl_exp(fpkm)
merge_df <- merge(gr_df, GL_exp[,c(1,7,20)], by.x='gene_id', by.y='ENSEMBL', all.x =T) 
sub_genes <- makeGRangesFromDataFrame(merge_df, keep.extra.columns=T)
### add 'chr' in subgenes
seqlevels(sub_genes) <-  paste('chr', seqlevels(sub_genes), sep='')
### make sub_genes_tss
sub_genes_tss <- promoters(sub_genes, upstream = 0, downstream = 1)
### Input Chip-seq
library(ChIPpeakAnno)
library(ChIPseeker)
ATAC <- readPeakFile('../CTRL_ATAC_S36_peaks.broadPeak', head=F)
K4 <- readPeakFile('../WT_K4_S24_macs2_peaks.narrowPeak', head=F)
K9 <- readPeakFile('../WT_K9_S25_macs2_peaks.broadPeak', head=F)
K27 <- readPeakFile('../WT_K27_S26_macs2_peaks.broadPeak', head=F)
H2AZ <- readPeakFile('../WTB-H2AZ_macs2_peaks.narrowPeak', head=F) 
### load WT methylation data
rrbs <- read.delim('~/MM_WT_RRBS.txt', header=T)
# convert it to GRanges
RRBS <- makeGRangesFromDataFrame(rrbs, keep.extra.columns = T,  seqnames.field = 'chr',
                                 strand.field="strand", start.field='start',
                                 end.field='end', starts.in.df.are.0based = F)
GL <- GRangesList(gr1=ATAC, gr2=K4, gr3=K9, gr4=K27, gr5=H2AZ, gr6=RRBS)
library(EnrichedHeatmap)
ma <- list()
for (i in c(1:5)){
  ma[[i]] <- normalizeToMatrix(GL[[i]], sub_genes_tss, value_column = "V5", smooth = T,
                                   extend = 5000, mean_mode = "w0", w = 50, background = 0)
}
### methylation data differs from Chip data
ma[[6]] <- normalizeToMatrix(RRBS, sub_genes_tss, value_column="methyl", smooth = T,
                                 extend = 5000, mean_mode = "absolute", w = 50,  background = 0)
### divided by 6 states
library(circlize)
col_fun = colorRamp2(c(0, 0.99), c("blue", "red"))
lgd = Legend(at = c('ES1', 'ES2', 'ES3', 'ES4', 'ES5', 'ES6'), title = "states", 
             type = "lines", legend_gp = gpar(col=1:6))
partition = sub_genes_tss$score
ht_list = Heatmap(partition, col=structure(1:6), 
                                           name = "state", show_row_names = FALSE, 
                                           width = unit(6, "mm")) +
  EnrichedHeatmap(ma[[1]], name="ATAC", column_title='ATAC', col=col_fun,
                  top_annotation=HeatmapAnnotation(lines=anno_enriched(gp=gpar(col=1:6)))) +
  EnrichedHeatmap(ma[[2]], name="H3K4me3", column_title='H3K4me3',col=col_fun,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp=gpar(col=1:6)))) +          
  EnrichedHeatmap(ma[[5]], name = "H2AZ", column_title = 'H2AZ',col=col_fun, 
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp=gpar(col=1:6)))) +          
  EnrichedHeatmap(ma[[3]], name="H3K9me3", column_title='H3K9me3', col=col_fun, 
                  top_annotation=HeatmapAnnotation(lines=anno_enriched(gp = gpar(col=1:6)))) +
  EnrichedHeatmap(ma[[4]], name = "H3K27me3", column_title = 'H3K27me3', col=col_fun, 
                  top_annotation=HeatmapAnnotation(lines = anno_enriched(gp=gpar(col=1:6)))) +
  EnrichedHeatmap(ma[[6]], name = "Methylation", column_title = 'Methylation', col=col_fun, 
                  top_annotation=HeatmapAnnotation(lines = anno_enriched(gp=gpar(col=c(1:6)))))
draw(ht_list, split = partition, annotation_legend_list = list(lgd))

