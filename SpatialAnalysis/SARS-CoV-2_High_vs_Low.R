source("./de_analysis.R")

#############
# Load data #
#############
wta.counts <- read.delim("data/Broad-COVID_WTA_Q3Norm_TargetCountMatrix.txt", row.names=1)

df_segments <- read.delim("data/Broad-COVID_WTA_SegmentProperties.txt", row.names=1)
rownames(df_segments) <- str_replace_all(rownames(df_segments), '-', '.')

df_tissue <- read.delim("data/annotation_file_wta.txt", row.names=1)
rownames(df_tissue) <- str_replace_all(rownames(df_tissue), '-', '.')
df_segments$tissue <- 'unknown'
df_segments[rownames(df_tissue), 'tissue'] <- df_tissue$Primary_Morph

# Extract donor info
df_segments$donor <- sapply(df_segments$scan.name, USE.NAMES=F, function(x) {
    l <- strsplit(x, '-')[[1]]
    return(l[1])
})

############################
# Assign SARS-CoV-2 Levels #
############################
df_segments$covid <- 'None'

df_covid.panck <- read.csv("WTA_SARS-CoV-2_PanCK.sig.csv", row.names=1, header=F)
rownames(df_covid.panck) <- str_replace_all(rownames(df_covid.panck), '-', '.')
colnames(df_covid.panck) <- c('sig_score')
df_covid.panck$covid <- sapply(df_covid.panck$sig_score, function(x) {
    if (x > 1.30) {
        return('High')
    } else if (x < 0.75) {
        return('Low')
    } else {
        return('Medium')
    }
})
df_segments[rownames(df_covid.panck),]$covid <- df_covid.panck$covid

df_covid.syto13 <- read.csv("WTA_SARS-CoV-2_Syto13.sig.csv", row.names=1, header=F)
rownames(df_covid.syto13) <- str_replace_all(rownames(df_covid.syto13), '-', '.')
colnames(df_covid.syto13) <- c('sig_score')
df_covid.syto13$covid <- sapply(df_covid.syto13$sig_score, function(x) {
    if (x > 1.30) {
        return('High')
    } else if (x < 0.75) {
        return('Low')
    } else {
        return('Medium')
    }
})
df_segments[rownames(df_covid.syto13),]$covid <- df_covid.syto13$covid

# Select only Alveolar AOIs
covid.segments <- subset(df_segments, segment %in% c("PanCK", "Syto13") & tissue=='Alveolar' & covid!='None',
                       select=c('slide.name', 'scan.name', 'roi', 'segment', 'covid', 'aoi', 'area', 'RawReads', 'NormFactorQ3'))

if (!dir.exists("Figure_4")) {
    dir.create("Figure_4")
}

if (!dir.exists("ED_Figure_8")) {
    dir.create("ED_Figure_8")
}

####################################
# Differential Expression Analysis #
####################################

## PanCK+
top.table11 <- analyze_by_contrasts(covid.segments, wta.counts, c('PanCK'), 'covid', c('Low', 'Medium', 'High'))
gen_de_csv(top.table11, "Figure_4/Figure_4e_DE")

pdf("Figure_4/Figure_4e.pdf", width=10, height=10)
EnhancedVolcano(
    top.table11,
    lab = rownames(top.table11),
    x = 'logFC',
    y = 'P.Value',
    selectLab = top_genes(top.table11),
    pCutoff = determine_p_cutoff(top.table11),
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
dev.off()

## PanCK-
top.table12 <- analyze_by_contrasts(covid.segments, wta.counts, c('Syto13'), 'covid', c('Low', 'Medium', 'High'))
gen_de_csv(top.table12, "ED_Figure_8/ED_Figure_8m_DE")

pdf("ED_Figure_8/ED_Figure_8m.pdf", width=10, height=10)
EnhancedVolcano(
    top.table12,
    lab = rownames(top.table12),
    x = 'logFC',
    y = 'P.Value',
    selectLab = top_genes(top.table12, logFC_thresh=0),
    pCutoff = determine_p_cutoff(top.table12),
    FCcutoff = 0,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
dev.off()
