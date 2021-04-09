source("./de_analysis.R")

#################
# Load CTA data #
#################
cta.counts <- read.delim("data/Broad-COVID_CTA_Q3Norm_TargetCountMatrix.txt", row.names=1)

df_segments <- read.delim("data/Broad-COVID_CTA_SegmentProperties.txt", row.names=1)
rownames(df_segments) <- str_replace_all(rownames(df_segments), '-', '.')

df_tissue <- read.delim("data/annotation_file_cta.txt", row.names=1)
rownames(df_tissue) <- str_replace_all(rownames(df_tissue), '-', '.')
df_segments$tissue <- 'unknown'
df_segments[rownames(df_tissue), 'tissue'] <- df_tissue$Primary_Morph

# Extract donor info
df_segments$donor <- sapply(df_segments$scan.name, USE.NAMES=F, function(x) {
    if (is.na(x)) {
        'unknown'
    } else {
        l <- strsplit(x, '_')[[1]]
        if (length(l) > 1) {
            l[1]
        } else {
            l <- strsplit(x, '-')[[1]]
            l[1]
        }
    }
})

# Split into Patient (S) and Control (C) groups
df_segments$group <- sapply(df_segments$donor, USE.NAMES=F, function(x) {
    if (is.na(x)) {
        'unknown'
    } else {
        substring(x, 1, 1)
    }
})

# Only use Alveolar AOIs
cta.segments <- subset(df_segments, donor!='unknown' & tissue=='Alveolar' & segment %in% c('PanCK', 'Syto13'),
                   select=c('slide.name', 'scan.name', 'roi', 'segment', 'group', 'tissue', 'donor', 'aoi', 'area', 'RawReads', 'NormFactorQ3'))

if (!dir.exists("ED_Figure_8")) {
    dir.create("ED_Figure_8")
}

####################################
# Differential Expression Analysis #
####################################

## PanCK+
top.table21 <- analyze_by_contrasts(cta.segments, cta.counts, c('PanCK'), 'group', c('C', 'S'))
gen_de_csv(top.table21, "ED_Figure_8/ED_Figure_8j_left_DE")

pdf("ED_Figure_8/ED_Figure_8j_left.pdf", width=10, height=10)
EnhancedVolcano(top.table21,
    lab = rownames(top.table21),
    x = 'logFC',
    y = 'P.Value',
    selectLab = top_genes(top.table21),
    pCutoff = determine_p_cutoff(top.table21),
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
dev.off()

## PanCK-
top.table22 <- analyze_by_contrasts(cta.segments, cta.counts, c('Syto13'), 'group', c('C', 'S'))
gen_de_csv(top.table22, "ED_Figure_8/ED_Figure_8k_left_DE")

pdf("ED_Figure_8/ED_Figure_8k_left.pdf", width=10, height=10)
EnhancedVolcano(top.table22,
    lab = rownames(top.table22),
    x = 'logFC',
    y = 'P.Value',
    selectLab = top_genes(top.table22, extra_genes=c('S', 'ORF1ab')),
    pCutoff = determine_p_cutoff(top.table22),
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
dev.off()

#######################################
# Gene Set Enrichment Analysis (GSEA) #
#######################################
gen_pathway_csv(top.table21, "ED_Figure_8/ED_Figure_8j_gsea", seed=0)
gen_pathway_csv(top.table22, "ED_Figure_8/ED_Figure_8k_gsea", seed=0)
