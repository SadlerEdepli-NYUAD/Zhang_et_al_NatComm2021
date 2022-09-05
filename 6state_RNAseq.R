library(genomation)
library(pheatmap)
library('biomaRt')
library(IRanges)
library(GenomicRanges)
library(ggplot2)

### seperate in-house 6 states
# add strand as '*' into column 4 before readBed
ES_file <- readBed('./Mouse_liver_6_segments.bed')
head(ES_file)
ES1 <- ES_file[grepl('E1', ES_file$score), ]
ES2 <- ES_file[ES_file$score == 'E2', ]
ES3 <- ES_file[ES_file$score == 'E3', ]
ES4 <- ES_file[ES_file$score == 'E4', ]
ES5 <- ES_file[ES_file$score == 'E5', ]
ES6 <- ES_file[ES_file$score == 'E6', ]
Genome_ES <- c(ES1, ES2, ES3, ES4, ES5, ES6)
# library("exomeCopy")
# ES1_sub <- subdivideGRanges(ES1, subsize=10)
# gl_sub <- GRangesList(ES1_sub, ES2_sub, ES3_sub, ES4_sub, ES5_sub, ES6_sub)
### pie chart summary
### compute genomic regions covered
### by base pair
slices <- c(sum(width(ES1)), sum(width(ES2)), sum(width(ES3)), 
            sum(width(ES4)), sum(width(ES5)), sum(width(ES6)))
per <- paste(round(100*slices/sum(slices), 1), '%', sep='')
lbls <- c("E1 Active TSS", "E2 Flanking TSS", 'E3 Transcription across genes',
          'E4 Heterochromatin', 'E5 Repetitive', 'E6 Repressive')
library(RColorBrewer)
setcols <- brewer.pal(6, "Set3")
pie(slices, labels = per, col = setcols) 
legend("topright", lbls, cex=0.6, fill=setcols) 
### Gene Annotation 
# readTranscriptFeatures(location, remove.unsual = TRUE, up.flank = 1000,
# down.flank = 1000, unique.prom = TRUE)
# set promoters as tss-500bp
gene.parts <- readTranscriptFeatures('mm10_refseq.bed', up.flank = 500, down.flank = 0)
gl <- GRangesList(ES1, ES2, ES3, ES4, ES5, ES6)
names(gl) <- lbls
anno_parts <- annotateWithGeneParts(gl, gene.parts, intersect.chr=T)
anno_genome <- annotateWithGeneParts(Genome_ES, gene.parts, intersect.chr=T)
anno_sub <- annotateWithGeneParts(gl_sub, gene.parts, intersect.chr=T)

library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
annotateWithGeneParts(genome, gene.parts, intersect.chr=T)
anno_df <- data.frame('Genome-wide'=c(9.25, 12.25, 22.91, 55.60),
                      'ES1'=c(13.09, 15.48, 46.56, 24.87),
                      'ES2'=c(34.04, 5.99, 23.39, 36.58),
                      'ES3'=c(4.34, 4.34, 25.65, 65.67),
                      'ES4'=c(7.08, 19.48, 16.23, 57.22),
                      'ES5'=c(0.40, 2.38, 21.04, 76.19),
                      'ES6'=c(13.98, 14.39, 20.43, 51.20))
rownames(anno_df) <- c('promoter', 'exon', 'intron', 'intergenic')
colnames(anno_df) <- c('Genome-wide', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6')
pal <- brewer.pal(4, "Set2")
barplot(as.matrix(anno_df), col=pal)
legend("topright", lbls, cex=1.0, fill=pal) 
### bp wised == proportion * sum(width(S)) 
c(13.09, 15.48, 46.56, 24.87) * sum(width(gl[[1]]))/100
=========================================================================================
### expression
wt_fpkm <- read.csv('./WT_FPKM_mean.csv', sep='\t', header=T)
### Sig_DEG_afterPH
library(xlsx)
sig_genes <- c()
for ( i in 1:8){ 
  dat <- read.xlsx('./Sig_DEG_afterPH_eachT.xlsx', sheetIndex = i)
  sig_genes <- c(sig_genes, as.character(dat$ensembl))
  sig_genes <- unique(sig_genes) # 2179 sig genes in 7 comparision
}

targets <- list('state1'=ES1, 'state2'=ES2, 'state3'=ES3, 
                'state4'=ES4, 'state5'=ES5, 'state6'=ES6)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library("org.Mm.eg.db")
library(ChIPseeker)
### annotate regions of each state
ES <- list(ES1, ES2, ES3, ES4, ES5, ES6)
Inform <- list()
### Annotation, tss into states
for (i in 1:6){
  Inform[[i]] <- as.data.frame(annotatePeak(ES[[i]], 
                                            TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                            annoDb = "org.Mm.eg.db"))
  Inform[[i]] <- subset(Inform[[i]], 
                        Inform[[i]]$distanceToTSS == 0 )
}
### remove duplicates, because distanceToTSS is for trannscript
Inform_df <- rbind(Inform[[1]], Inform[[2]], Inform[[3]],
                   Inform[[4]], Inform[[5]], Inform[[6]])
library(dplyr)
Inform_df_filter <- Inform_df %>% 
  group_by(ENSEMBL) %>% 
  filter(transcriptId==min(transcriptId))
### rewrite the Inform list without transcript duplicates
Inform[[1]] <- Inform_df_filter %>% filter(score=='E1')
Inform[[2]] <- Inform_df_filter %>% filter(score=='E2')
Inform[[3]] <- Inform_df_filter %>% filter(score=='E3')
Inform[[4]] <- Inform_df_filter %>% filter(score=='E4')
Inform[[5]] <- Inform_df_filter %>% filter(score=='E5')
Inform[[6]] <- Inform_df_filter %>% filter(score=='E6')
 
head(Inform[[2]])
### build a FPKM list EXP and a gene state list GS with ES state
EXP <- list()
GS <- list()
for (i in 1:6){
  ### limit annotation within TSS +- 100 bp
  # inform_in <- Inform[[i]][which(Inform[[i]]$distanceToTSS == 0), ]
  # inform_2kb <- Inform[[i]][which(Inform[[i]]$distanceToTSS <= 100 &
  #                                  Inform[[i]]$distanceToTSS >= -100), ]
                                    #  Inform[[i]]$distanceToTSS >= -1000), ]
  # inform_2kb for defining the ES of a gene
  # get gene id, 17: ennsembl, 18: symbol
  # id <- inform_1kb[ ,c(17,18)]
  # 6: is ES 
  id <- Inform[[i]][ ,c(6,17,18)]
  # drop duplicates
  id <- id[!duplicated(id$ENSEMBL),]
  GS[[i]] <- id
  # get fpkm from whole data
  fpkm <- wt_fpkm[which(rownames(wt_fpkm) %in% id$ENSEMBL), ] 
  # filter expressed genes: fpkm>10, across all time points, 
  # which means genes at least expressed at some point during LR
  # fpkm <- fpkm[rowSums(fpkm)>10, ]
  EXP[[i]] <- fpkm 
}
====================================================================================
### get .bed files with ensembl id by using BioMart
library('biomaRt')
listMarts()    # to see which database options are present
ensembl <- useMart('ensembl')
mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))

S1_Gene_bed <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position',
                  'strand', 'ensembl_gene_id', 'mgi_symbol'),
      filters='ensembl_gene_id',values=GS[[1]]$ENSEMBL, mart=mart_mm) 
### add 'chr' to chromosome name, change strand
S1_Gene_bed$chromosome_name <- paste('chr', S1_Gene_bed$chromosome_name, sep = '')
S1_Gene_bed$strand <-  with(S1_Gene_bed, 
                            ifelse(strand=='1', '+', '-'))
write.table(S1_Gene_bed, 'S1_gene.txt', row.names=F, quote=F, 
            col.names=F, sep='\t')
===============================================================================================
### Gene expression in 6 states
### make the matrix for plotting
df <- data.frame(fpkm = c(EXP[[1]]$baselineWT_1, EXP[[2]]$baselineWT_1,
                           EXP[[3]]$baselineWT_1, EXP[[4]]$baselineWT_1,
                           EXP[[5]]$baselineWT_1, EXP[[6]]$baselineWT_1),
                  state = c(rep('ES1', nrow(EXP[[1]])), rep('ES2', nrow(EXP[[2]])),
                            rep('ES3', nrow(EXP[[3]])), rep('ES4', nrow(EXP[[4]])),
                            rep('ES5', nrow(EXP[[5]])), rep('ES6', nrow(EXP[[6]]))),
                  ensembl = c(rownames(EXP[[1]]), rownames(EXP[[2]]),
                              rownames(EXP[[3]]), rownames(EXP[[4]]),
                              rownames(EXP[[5]]), rownames(EXP[[6]])))
### set fpkm > 200 == 200 for plotting
df_p200 <- df
df_p200$fpkm <- ifelse(df_p200$fpkm > 200, 200, df_p200$fpkm)
p <- ggplot(df_p200, aes(state, fpkm))+ geom_jitter(alpha = I(1 / 2)) + ylim(0, 200)
p
### ============input liver genes
lvr_gene <- read.delim('../Curated_LiverGenes.txt')
lvr_enriched_gene <- read.delim('../Liver_TissueEnriched.tsv')
### convert hg lvr gene to mm ensembl ID
ensembl <- useMart('ensembl')  # using ensembl database data
mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
mart_hg <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))

mm_lvr_gene <- getLDS(filters= "ensembl_gene_id", attributes= "ensembl_gene_id", 
                      mart =mart_hg, values=lvr_enriched_gene$Ensembl, 
                      attributesL = "ensembl_gene_id", martL =  mart_mm)

df_p200$lvr <- df_p200$ensembl %in% unique(mm_lvr_gene$Gene.stable.ID.1)
p <- ggplot(df_p200, aes(state, fpkm, col=sig)) + geom_jitter(alpha = I(1 / 2)) + 
  ylim(0, 200) + scale_color_manual(values=c("#999999","red")) + theme_classic()
p_lvr <- ggplot(df_p200, aes(state, fpkm, col=lvr)) + geom_jitter(alpha = I(1 / 2)) + 
  ylim(0, 200) + scale_color_manual(values=c("#999999","red")) + theme_classic()
p_lvr
### bar plot to show expressed and non-expressed genes
### baseline fpkm >= 1 defined as expressed
dat_Num_Exp <- data.frame('ES1'= as.vector(table(df_p200[df_p200$state=='ES1', 1] >= 1)), 
                          'ES2'= as.vector(table(df_p200[df_p200$state=='ES2', 1] >= 1)),
                          'ES3'= as.vector(table(df_p200[df_p200$state=='ES3', 1] >= 1)),
                          'ES4'= as.vector(table(df_p200[df_p200$state=='ES4', 1] >= 1)),
                          'ES5'= as.vector(table(df_p200[df_p200$state=='ES5', 1] >= 1)),
                          'ES6'= as.vector(table(df_p200[df_p200$state=='ES6', 1] >= 1)))
rownames(dat_Num_Exp) <- c('Silenced', 'Expressed')
barplot(as.matrix(dat_Num_Exp), ylim = c(1000,8000))
### ======= bivalent genes=============================================
#k27_genes <- read.table('../Regeneration/Sub_Gene_Set/K27_positive_genes.bed')
#library('biomaRt')
#ensembl <- useMart('ensembl')  # using ensembl database data
#mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))
#k27_genes_ensembl <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
#                          filters='entrezgene_id',values=k27_genes$V4, mart=mart_mm) 
#head(k27_genes_ensembl)
#df_p200$k27 <- df_p200$ensembl %in% unique(k27_genes_ensembl$ensembl_gene_id)
#p_k27 <- ggplot(df_p200, aes(state, fpkm, col=k27)) + geom_jitter(alpha = I(1 / 2)) + 
#  ylim(0, 200) + scale_color_manual(values=c("#999999","red")) + theme_classic()
#p_k27
K27_WT_2nd <- readPeakFile('../WT_K27_S26_macs2_peaks.broadPeak', head=F)

k27_anno <- as.data.frame(annotatePeak(K27_WT_2nd, 
                                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                          annoDb = "org.Mm.eg.db"))
k27_anno <- subset(k27_anno, k27_anno$distanceToTSS == 0 )
k27_anno <- k27_anno %>% 
  group_by(ENSEMBL) %>% 
  filter(transcriptId==min(transcriptId))
# colnames(k27_anno): 21: ensembl; 22: gene symbol
k27_id <- k27_anno[ ,c(21, 22)]


df_p200$k27 <- df_p200$ensembl %in% k27_id$ENSEMBL
p_k27 <- ggplot(df_p200, aes(state, fpkm, col=k27)) + geom_jitter(alpha = I(1 / 2)) + 
  ylim(0, 200) + scale_color_manual(values=c("#999999","blue")) + theme_classic()
p_k27
# boxplot state1 & k27+ state1
df1 <- df[df$state=='ES1', ]
df1$k27 <- df1$ensembl %in% k27_id$ENSEMBL
p_bp <- ggplot(df1, aes(x=k27, y=log(fpkm+1))) + 
  geom_violin()
p_bp
================================================================================================
### df is with fpkm, state and ensembl id
head(df)
### make df as a GRagnen window by integrating with anotation ensembl_gene
library(ensembldb)
library(EnsDb.Mmusculus.v79)
names(listTables(EnsDb.Mmusculus.v79))
edb <- EnsDb.Mmusculus.v79  # abbreviate
ensembl_genes <- genes(edb)
ensembl_genes
sub_genes <-  subset(ensembl_genes, gene_id %in% df$ensembl)

#====6 States GO Analysis=================================================================================
library(clusterProfiler)
library(org.Mm.eg.db) ### divide into expressed and silenced genes 
### barplot threshold 1
E1_gene <- filter(df, state=='ES1' & fpkm > 5)
E1_gene <- filter(df, state=='ES1' & fpkm <= 5)
E2_gene <- filter(df, state=='ES2' & fpkm > 5)
E2_gene <- filter(df, state=='ES2' & fpkm <= 5)
E3_gene <- filter(df, state=='ES3' & fpkm > 5)
E3_gene <- filter(df, state=='ES3' & fpkm <= 5)
E4_gene <- filter(df, state=='ES4' & fpkm > 5)
E4_gene <- filter(df, state=='ES4' & fpkm <= 5)
E5_gene <- filter(df, state=='ES5')
E6_gene <- filter(df, state=='ES6')

E1_df <- bitr(E1_gene$ensembl, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
E2_df <- bitr(E2_gene$ensembl, fromType = "ENSEMBL",
              toType = c("ENTREZID", "SYMBOL"),
              OrgDb = org.Mm.eg.db)
E3_df <- bitr(E3_gene$ensembl, fromType = "ENSEMBL",
              toType = c("ENTREZID", "SYMBOL"),
              OrgDb = org.Mm.eg.db)
E4_df <- bitr(E4_gene$ensembl, fromType = "ENSEMBL",
              toType = c("ENTREZID", "SYMBOL"),
              OrgDb = org.Mm.eg.db)
E5_df <- bitr(E5_gene$ensembl, fromType = "ENSEMBL",
              toType = c("ENTREZID", "SYMBOL"),
              OrgDb = org.Mm.eg.db)
E6_df <- bitr(E6_gene$ensembl, fromType = "ENSEMBL",
              toType = c("ENTREZID", "SYMBOL"),
              OrgDb = org.Mm.eg.db)
Clusterlist <- list(E1_df$ENTREZID, E2_df$ENTREZID)
names(Clusterlist) <- c('S1', 'S2')
str(Clusterlist)
CompareGO_BP <- compareCluster(Clusterlist, fun="enrichGO", 
                               pvalueCutoff=0.01, pAdjustMethod="BH", 
                               OrgDb=org.Mm.eg.db, ont="BP", readable=T)

CompareGO_MF <- compareCluster(Clusterlist, fun="enrichGO", 
                               pvalueCutoff=0.01, pAdjustMethod="BH", 
                               OrgDb=org.Mm.eg.db, ont="MF", readable=T)


dotplot(CompareGO_BP, showCategory=10, title="GO - Biological Process")
dotplot(CompareGO_MF, showCategory=10, title="GO - Molecular Function")

E4_ego <- clusterProfiler::enrichGO(gene     = E4_df$ENTREZID,
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05, 
                                 readable      = TRUE)
barplot(E2_ego, showCategory=20)
write.table(as.data.frame(E4_ego)[1:50, ], 'E4_Silenced_GO.txt', quote = F, sep='\t')
#========================================================================================================
### generate gene list for the heatmap
### cell cycle gene list
Cell_cyc <- read.table('../Regeneration/K27_TC/Cellcycle_gene.bed')
### Input S1 and S2 silenced genes
S1_neg <- read.table('../Regeneration/Figures_progress/Fig3/GO_results/E1_Silenced_GO.txt', sep='\t' )
### select certain pathways
S1_neg <- subset(S1_neg, Description=='epithelial cell proliferation' | 
                   Description=='regulation of epithelial cell proliferation')
### narrow downn to the col
S1_neg_ID <- S1_neg$geneID
### split by '/' 
S1_neg_ID <- strsplit(as.character(S1_neg_ID), '/')
### unlist, then unique
S1_neg_ID <- unique(unlist(S1_neg_ID))

### Input S1 and S2 silenced genes
S2_neg <- read.table('../Regeneration/Figures_progress/Fig3/GO_results/E2_Silenced_GO.txt', sep='\t' )
S2_neg <- subset(S2_neg, Description=='DNA repair' | 
                   Description=='double-strand break repair' |
                   Description=='mitotic cell cycle phase transition' |
                   Description=='cell cycle phase transition' |
                   Description=='DNA recombination' |
                   Description=='double-strand break repair via homologous recombination'|
                   Description=='recombinational repair' |
                   Description=='chromosome segregation' |
                   Description=='nuclear division' |
                   Description=='regulation of DNA metabolic process' |
                   Description=='regulation of chromosome organization' |
                   Description=='DNA replication' |
                   Description=='regulation of mitotic cell cycle')
### narrow downn to the col
S2_neg_ID <- S2_neg$geneID
### split by '/' 
S2_neg_ID <- strsplit(as.character(S2_neg_ID), '/')
### unlist, then unique
S2_neg_ID <- unique(unlist(S2_neg_ID))
### combine S1 and S2
neg_ID <- c(S1_neg_ID, S2_neg_ID)

library('biomaRt')
listMarts()    # to see which database options are present
ensembl <- useMart('ensembl')
mart_mm <- useDataset('mmusculus_gene_ensembl', useMart('ensembl'))

S12_neg_Gene_bed <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position',
                                  'strand', 'ensembl_gene_id', 'mgi_symbol'),
                     filters='mgi_symbol',values=neg_ID, mart=mart_mm) 
### add 'chr' to chromosome name, change strand
S12_neg_Gene_bed$chromosome_name <- paste('chr', S12_neg_Gene_bed$chromosome_name, 
                                          sep = '')
S12_neg_Gene_bed$strand <-  with(S12_neg_Gene_bed, 
                            ifelse(strand=='1', '+', '-'))
### output .txt file
write.table(S12_neg_Gene_bed, 'S12_neg_Gene.txt', row.names=F, quote=F, 
            col.names=F, sep='\t')
### output .bed file
df12 <- data.frame(seqnames = S12_neg_Gene_bed$chromosome_name,
                  starts = S12_neg_Gene_bed$start_position-1,
                  ends = S12_neg_Gene_bed$end_position, 
                  scores = S12_neg_Gene_bed$ensembl_gene_id,
                  names = S12_neg_Gene_bed$mgi_symbol,
                  strands = S12_neg_Gene_bed$strand)
write.table(df12, file="S12_neg_Gene.bed", quote=F, sep="\t", row.names=F, col.names=F)

