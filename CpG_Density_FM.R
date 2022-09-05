# find all genome CpGs in Mouse (all the chromosome including mitocondrial)
########################################################################################################################################
library(BSgenome.Mmusculus.UCSC.mm10)  
chrs <- names(Mmusculus)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 2))))
#########################################################################################################

# import and manipulate Repetitive Elements locations
#########################################################################################################
All_TE <- read.delim("./Input_Data/Mouse_mm10_rmsk.gz")
head(All_TE)

# reodering columns
All_TE = All_TE[c("genoName", "genoStart", "genoEnd", "strand",
                  "repName", "repClass", "repFamily", "swScore",
                  "milliDiv", "milliDel", "milliIns", "X.bin",
                  "repStart", "repEnd", "repLeft", "genoLeft", "id")]

# remove unused columns
All_TE <- subset(All_TE, select = -c(milliDiv, milliDel, milliIns,
                                     X.bin, repStart, repEnd,
                                     repLeft, genoLeft, id))

# renaming columns for conversion in Genomic Ranges
colnames(All_TE) <- c("chr", "start", "end", "strand", "repName",
                      "repClass", "repFamily", "swScore")
#########################################################################################################

# Import H3K9me3 broadPeaks file as GRanges
########################################################################################################################################
library(genomation)
bed.file = "./Input_Data/WT_K9_S25_macs2_peaks.broadPeak"
H3K9me3_peaks_MouseLiver_Control_GR = readBed(bed.file)
colnames(mcols(H3K9me3_peaks_MouseLiver_Control_GR)) <- c("score", "name", "signalValue", "pValue", "qValue")
H3K9me3_peaks_MouseLiver_Control_GR
########################################################################################################################################

# Import States bed files as GRanges
########################################################################################################################################
library(genomation)
bed.file = "./Input_Data/State1.bed"
State_1_GR = readBed(bed.file)
colnames(mcols(State_1_GR)) <- c("State_ID", "strand")
bed.file = "./Input_Data/State2.bed"
State_2_GR = readBed(bed.file)
colnames(mcols(State_2_GR)) <- c("State_ID", "strand")
bed.file = "./Input_Data/State3.bed"
State_3_GR = readBed(bed.file)
colnames(mcols(State_3_GR)) <- c("State_ID", "strand")
bed.file = "./Input_Data/State4.bed"
State_4_GR = readBed(bed.file)
colnames(mcols(State_4_GR)) <- c("State_ID", "strand")
bed.file = "./Input_Data/State5.bed"
State_5_GR = readBed(bed.file)
colnames(mcols(State_5_GR)) <- c("State_ID", "strand")
bed.file = "./Input_Data/State6.bed"
State_6_GR = readBed(bed.file)
colnames(mcols(State_6_GR)) <- c("State_ID", "strand")
#########################################################################################################

# importing RRBS CpGs
#########################################################################################################
CpGs_All <- read.delim("./Input_Data/CpGs_All.txt")
# converting in GRanges
GR_CpGs_All <- as(CpGs_All, "GRanges")
#########################################################################################################

# subsetting data.frame
#########################################################################################################
# all Transposons
All_TE_Transposons <- subset(All_TE, repClass==c("LTR") | repClass==("SINE") | repClass==c("LINE") | repClass==("DNA"))
# all Non-Transposons
All_TE_Others <- subset(All_TE, repClass!=c("LTR") & repClass!=("SINE") & repClass!=c("LINE") & repClass!=("DNA"))
# all Transposons dived by class (LINE, SINE, LTR and DNA)
All_TE_LTR <- subset(All_TE, repClass==("LTR"))
All_TE_SINE <- subset(All_TE, repClass==("SINE"))
All_TE_LINE <- subset(All_TE, repClass==("LINE"))
All_TE_DNA <- subset(All_TE, repClass==("DNA"))
# converting in GRanges
#########################################################################################################
All_TE_GR <- as(All_TE, "GRanges") # all Repetitive Elements
All_TE_Transposons_GR <- as(All_TE_Transposons, "GRanges") # all Transposons
All_TE_Others_GR <- as(All_TE_Others, "GRanges") # all Non-Transposons
All_TE_LTR_GR <- as(All_TE_LTR, "GRanges")
All_TE_SINE_GR <- as(All_TE_SINE, "GRanges")
All_TE_LINE_GR <- as(All_TE_LINE, "GRanges")
All_TE_DNA_GR <- as(All_TE_DNA, "GRanges")
#########################################################################################################

# subset By Overlaps GRanges
#########################################################################################################
CpGs_State_1_GR <- subsetByOverlaps(GR_CpGs_All, State_1_GR)
CpGs_State_2_GR <- subsetByOverlaps(GR_CpGs_All, State_2_GR)
CpGs_State_3_GR <- subsetByOverlaps(GR_CpGs_All, State_3_GR)
CpGs_State_4_GR <- subsetByOverlaps(GR_CpGs_All, State_4_GR)
CpGs_State_5_GR <- subsetByOverlaps(GR_CpGs_All, State_5_GR)
CpGs_State_6_GR <- subsetByOverlaps(GR_CpGs_All, State_6_GR)

K9_All_TE_GR <- subsetByOverlaps(All_TE_GR, H3K9me3_peaks_MouseLiver_Control_GR)
noK9_All_TE_GR <- subsetByOverlaps(All_TE_GR, H3K9me3_peaks_MouseLiver_Control_GR, invert=TRUE)
K9_All_TE_Transposons_GR <- subsetByOverlaps(All_TE_Transposons_GR, H3K9me3_peaks_MouseLiver_Control_GR)
noK9_All_TE_Transposons_GR <- subsetByOverlaps(All_TE_Transposons_GR, H3K9me3_peaks_MouseLiver_Control_GR, invert=TRUE)
K9_All_TE_Others_GR <- subsetByOverlaps(All_TE_Others_GR, H3K9me3_peaks_MouseLiver_Control_GR)
noK9_All_TE_Others_GR <- subsetByOverlaps(All_TE_Others_GR, H3K9me3_peaks_MouseLiver_Control_GR, invert=TRUE)

CpGs_All_TE_Transposons_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_GR)
CpGs_All_TE_Others_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Others_GR)
CpGs_All_TE_LTR_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_LTR_GR)
CpGs_All_TE_SINE_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_SINE_GR)
CpGs_All_TE_LINE_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_LINE_GR)
CpGs_All_TE_DNA_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_DNA_GR)

All_TE_Transposons_State_1_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_1_GR)
All_TE_Transposons_State_2_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_2_GR)
All_TE_Transposons_State_3_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_3_GR)
All_TE_Transposons_State_4_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_4_GR)
All_TE_Transposons_State_5_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_5_GR)
All_TE_Transposons_State_6_GR <- subsetByOverlaps(All_TE_Transposons_GR, State_6_GR)
CpGs_All_TE_Transposons_State_1_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_1_GR)
CpGs_All_TE_Transposons_State_2_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_2_GR)
CpGs_All_TE_Transposons_State_3_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_3_GR)
CpGs_All_TE_Transposons_State_4_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_4_GR)
CpGs_All_TE_Transposons_State_5_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_5_GR)
CpGs_All_TE_Transposons_State_6_GR <- subsetByOverlaps(GR_CpGs_All, All_TE_Transposons_State_6_GR)

K9_All_TE_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, K9_All_TE_GR)
noK9_All_TE_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, noK9_All_TE_GR)
K9_All_TE_Transposons_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, K9_All_TE_Transposons_GR)
noK9_All_TE_Transposons_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, noK9_All_TE_Transposons_GR)
K9_All_TE_Others_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, K9_All_TE_Others_GR)
noK9_All_TE_Others_CpGs_All_GR <- subsetByOverlaps(GR_CpGs_All, noK9_All_TE_Others_GR)
#########################################################################################################

#########################################################################################################
# All CpGs in States perc.meth boxplot Fig.2B
#########################################################################################################
# use lengths because is a dataframe already otherwise is length for GRanges
sq <- seq(max(length(CpGs_State_1_GR),
              length(CpGs_State_2_GR),
              length(CpGs_State_3_GR),
              length(CpGs_State_4_GR),
              length(CpGs_State_5_GR),
              length(CpGs_State_6_GR)
))

# this is forcing to keep the same lenght for all columns and put NAs
df_States_perc.meth_MouseLiver_Control <- data.frame(CpGs_State_1_GR$MouseLiver_Control[sq],
                                                     CpGs_State_2_GR$MouseLiver_Control[sq],
                                                     CpGs_State_3_GR$MouseLiver_Control[sq],
                                                     CpGs_State_4_GR$MouseLiver_Control[sq],
                                                     CpGs_State_5_GR$MouseLiver_Control[sq],
                                                     CpGs_State_6_GR$MouseLiver_Control[sq]
)

colnames(df_States_perc.meth_MouseLiver_Control) <- c("State_1",
                                                      "State_2",
                                                      "State_3",
                                                      "State_4",
                                                      "State_5",
                                                      "State_6"
)

library("reshape2")
df_States_perc.meth_MouseLiver_Control_melt <- melt(df_States_perc.meth_MouseLiver_Control)

library("ggplot2")
ggplot(df_States_perc.meth_MouseLiver_Control_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "perc.meth_MouseLiver_Control", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#F4C700", "#2D9037", "#E66654", "#B1529A", "#92CEC2", "#2D4596"))
#########################################################################################################

#########################################################################################################
# All Transposons - H3K9me3 positive or negative - Whole Genome CpGs density metaplot Fig.4F
#########################################################################################################
library(genomation)
ScoreBin_CpG_WG_K9_All_TE_Transposons = ScoreMatrixBin(target=cpgr,
                                                       windows=K9_All_TE_Transposons_GR,
                                                       bin.num=40, strand.aware=TRUE)

ScoreBin_CpG_WG_noK9_All_TE_Transposons = ScoreMatrixBin(target=cpgr,
                                                         windows=noK9_All_TE_Transposons_GR,
                                                         bin.num=40, strand.aware=TRUE)

# create list for merged plot
########################################################################
ScoreBin_CpG_WG_K9_noK9_All_Transposons_list <- list(ScoreBin_CpG_WG_K9_All_TE_Transposons, ScoreBin_CpG_WG_noK9_All_TE_Transposons)
ScoreBin_CpG_WG_K9_noK9_All_Transposons = ScoreMatrixList(target=ScoreBin_CpG_WG_K9_noK9_All_Transposons_list,
                                                          bin.num=40,
                                                          strand.aware=TRUE)

plotMeta(ScoreBin_CpG_WG_K9_noK9_All_Transposons,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.041), 
         profile.names=c("K9_All_Transposons", "noK9_All_Transposons"),
         winsorize=c(0,99),line.col=c("#127AB8", "darkgray"), lwd=2, lty = 1,
         ylab = "Avg. Score on CpGs density", xlab = "Window bins")
#########################################################################################################

# All Non Transposons - H3K9me3 positive or negative - Whole Genome CpGs density metaplot Fig.4F
#########################################################################################################
library(genomation)
ScoreBin_CpG_WG_K9_All_TE_Others = ScoreMatrixBin(target=cpgr,
                                                  windows=K9_All_TE_Others_GR,
                                                  bin.num=40, strand.aware=TRUE)

ScoreBin_CpG_WG_noK9_All_TE_Others = ScoreMatrixBin(target=cpgr,
                                                    windows=noK9_All_TE_Others_GR,
                                                    bin.num=40, strand.aware=TRUE)

# create list for merged plot
########################################################################
ScoreBin_CpG_WG_K9_noK9_All_TE_Others_list <- list(ScoreBin_CpG_WG_K9_All_TE_Others, ScoreBin_CpG_WG_noK9_All_TE_Others)
ScoreBin_CpG_WG_K9_noK9_All_TE_Others = ScoreMatrixList(target=ScoreBin_CpG_WG_K9_noK9_All_TE_Others_list,
                                                        bin.num=40,
                                                        strand.aware=TRUE)

plotMeta(ScoreBin_CpG_WG_K9_noK9_All_TE_Others,
         smoothfun=function(x) stats::lowess(x, f = 1/20),
         ylim = c(0, 0.041), 
         profile.names=c("K9_All_Non_Transposons", "noK9_All_Non_Transposons"),
         winsorize=c(0,99),line.col=c("#127AB8", "darkgray"), lwd=2, lty = 2,
         ylab = "Avg. Score on CpGs density", xlab = "Window bins")
#########################################################################################################

#########################################################################################################
# All Transposons and Non-Transposons - RRBS CpGs perc.meth boxplot Suppl Fig.3B
#########################################################################################################
# All MouseLiver_Control - All TE divided for repFamily
sq <- seq(max(length(CpGs_All_TE_DNA_GR),
              length(CpGs_All_TE_LTR_GR),
              length(CpGs_All_TE_LINE_GR),
              length(CpGs_All_TE_SINE_GR),
              length(CpGs_All_TE_Others_GR)))
# this is forcing to keep the same lenght for all columns and put NAs
df_perc.meth_All_TE_groups_MouseLiver_Control <- data.frame(CpGs_All_TE_DNA_GR$MouseLiver_Control[sq],
                                                            CpGs_All_TE_LTR_GR$MouseLiver_Control[sq],
                                                            CpGs_All_TE_LINE_GR$MouseLiver_Control[sq],
                                                            CpGs_All_TE_SINE_GR$MouseLiver_Control[sq],
                                                            CpGs_All_TE_Others_GR$MouseLiver_Control[sq])
head(df_perc.meth_All_TE_groups_MouseLiver_Control)
tail(df_perc.meth_All_TE_groups_MouseLiver_Control)
colnames(df_perc.meth_All_TE_groups_MouseLiver_Control) <- c("DNA",
                                                             "LTR",
                                                             "LINE",
                                                             "SINE",
                                                             "Non-TEs")

library("reshape2")
df_perc.meth_All_TE_groups_MouseLiver_Control_melt <- melt(df_perc.meth_All_TE_groups_MouseLiver_Control)

library("ggplot2")
ggplot(df_perc.meth_All_TE_groups_MouseLiver_Control_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "All TE families", y = "Percentage of CpGs methylation values") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#5B348C", "#D3A8CF", "#B54592", "#5D2262", "#B5B5B4"))
#########################################################################################################

#########################################################################################################
# All CpGs in Transposons divided by States - RRBS CpGs perc.meth boxplot Suppl Fig.3C and D
#########################################################################################################
# use lengths because is a dataframe already otherwise is length for GRanges
sq <- seq(max(length(CpGs_All_TE_Transposons_State_1_GR),
              length(CpGs_All_TE_Transposons_State_2_GR),
              length(CpGs_All_TE_Transposons_State_3_GR),
              length(CpGs_All_TE_Transposons_State_4_GR),
              length(CpGs_All_TE_Transposons_State_5_GR),
              length(CpGs_All_TE_Transposons_State_6_GR)
))

# this is forcing to keep the same lenght for all columns and put NAs
df_All_TE_Transposons_States_perc.meth_MouseLiver_Control <- data.frame(CpGs_All_TE_Transposons_State_1_GR$MouseLiver_Control[sq],
                                                                        CpGs_All_TE_Transposons_State_2_GR$MouseLiver_Control[sq],
                                                                        CpGs_All_TE_Transposons_State_3_GR$MouseLiver_Control[sq],
                                                                        CpGs_All_TE_Transposons_State_4_GR$MouseLiver_Control[sq],
                                                                        CpGs_All_TE_Transposons_State_5_GR$MouseLiver_Control[sq],
                                                                        CpGs_All_TE_Transposons_State_6_GR$MouseLiver_Control[sq]
)

colnames(df_All_TE_Transposons_States_perc.meth_MouseLiver_Control) <- c("All_Transposons_State_1",
                                                                         "All_Transposons_State_2",
                                                                         "All_Transposons_State_3",
                                                                         "All_Transposons_State_4",
                                                                         "All_Transposons_State_5",
                                                                         "All_Transposons_State_6"
)

library("reshape2")
df_All_TE_Transposons_States_perc.meth_MouseLiver_Control_melt <- melt(df_All_TE_Transposons_States_perc.meth_MouseLiver_Control)

library("ggplot2")
ggplot(df_All_TE_Transposons_States_perc.meth_MouseLiver_Control_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "perc.meth_MouseLiver_Control", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#F4C700", "#2D9037", "#E66654", "#B1529A", "#92CEC2", "#2D4596"))

ggplot(df_All_TE_Transposons_States_perc.meth_MouseLiver_Control_melt, aes(value, after_stat(count))) +
  geom_density(aes(color = variable, fill=variable), size= 1, alpha=0.25) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  labs(y= "Count (Density*Number) of CpGs per each sample", x= "% of methylation") +
  theme_classic() + scale_fill_manual(values = c("#F4C700", "#2D9037", "#E66654", "#B1529A", "#92CEC2", "#2D4596")) +
  scale_color_manual(values=c("#F4C700", "#2D9037", "#E66654", "#B1529A", "#92CEC2", "#2D4596"))
#########################################################################################################

#########################################################################################################
# All CpGs in Repetive Elements (All-REs, TEs and Non-TEs) - H3K9me3 positive or negative - 
# RRBS CpGs perc.meth boxplot Suppl Fig.3E
#########################################################################################################
# use lengths because is a dataframe already otherwise is length for GRanges
sq <- seq(max(length(K9_All_TE_CpGs_All_GR),
              length(noK9_All_TE_CpGs_All_GR),
              length(K9_All_TE_Transposons_CpGs_All_GR),
              length(noK9_All_TE_Transposons_CpGs_All_GR),
              length(K9_All_TE_Others_CpGs_All_GR),
              length(noK9_All_TE_Others_CpGs_All_GR)
))

# this is forcing to keep the same lenght for all columns and put NAs
df_K9_All_TE_groups_perc.meth_MouseLiver_Control <- data.frame(K9_All_TE_CpGs_All_GR$MouseLiver_Control[sq],
                                                               noK9_All_TE_CpGs_All_GR$MouseLiver_Control[sq],
                                                               K9_All_TE_Transposons_CpGs_All_GR$MouseLiver_Control[sq],
                                                               noK9_All_TE_Transposons_CpGs_All_GR$MouseLiver_Control[sq],
                                                               K9_All_TE_Others_CpGs_All_GR$MouseLiver_Control[sq],
                                                               noK9_All_TE_Others_CpGs_All_GR$MouseLiver_Control[sq]
)

colnames(df_K9_All_TE_groups_perc.meth_MouseLiver_Control) <- c("All_RepetitiveElements_K9",
                                                                "All_RepetitiveElements_noK9",
                                                                "All_Transposons_K9",
                                                                "All_Transposons_noK9",
                                                                "All_notTransposons_K9",
                                                                "All_notTransposons_noK9"
)

library("reshape2")
df_K9_All_TE_groups_perc.meth_MouseLiver_Control_melt <- melt(df_K9_All_TE_groups_perc.meth_MouseLiver_Control)

library("ggplot2")
ggplot(df_K9_All_TE_groups_perc.meth_MouseLiver_Control_melt, aes(x = variable, y = value, fill=variable)) + 
  geom_boxplot(width=0.5, notch = FALSE) + scale_y_continuous(breaks=seq(0,100,25)) + 
  labs(x= "RE_groups", y = "% of methylation") + 
  theme_classic() + theme(axis.text.x = element_text(angle=30, hjust = 1)) + 
  scale_fill_manual(values = c("#127AB8", "darkgray", "#127AB8", "darkgray", "#127AB8", "darkgray"))
#########################################################################################################


