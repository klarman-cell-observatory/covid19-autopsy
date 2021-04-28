rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)
library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(genefilter)
library(GSVA)
library(GSEABase)
library(pheatmap)

# ===== Functions =====
GetMidPQuantileResiduals <- function(fit){
  # fit: a GLM fit from edgeR [DGEGLM]
  y <- fit$counts
  mu <- fit$fitted.values
  phi <- fit$dispersion
  
  res <- zscoreNBinom(y,mu,size=1/phi)
  return(res)
}

# ==== Lung Data =====
# sample annotation
lung_annot = read.table(
  file = "/data/work/Projects/BrScRNAseq/Plots/BroadPaper/CODE/data/LungDSP_sample_annotation.tsv",
  header = T, sep = "\t"
)

# normalized counts
lung_q3 = read.table(
  file = "/data/work/Projects/BrScRNAseq/Plots/BroadPaper/CODE/data/LungDSP_Q3Norm_counts.tsv",
  header = T, sep = "\t", row.names = 1
)

# ==== Normalization =====
lungDSP_counts <- DGEList(counts= lung_q3)
keep.exprs <- filterByExpr(lungDSP_counts, group=lung_annot$subtissue)
lungDSP_counts <- lungDSP_counts[keep.exprs,, keep.lib.sizes=FALSE]
# normalize counts (TMM)
lungDSP_counts <- calcNormFactors(object=lungDSP_counts,method="TMM")

lungDSP_cpm <- cpm(lungDSP_counts,log=T)

# ==== Adjust Expression =====
edgerres_model <- model.matrix(~ lung_annot$Patient)

# estimate dispersion
resid_dge = estimateGLMCommonDisp(y=lungDSP_counts,design = edgerres_model, verbose=TRUE)
resid_dge = estimateGLMTrendedDisp(y=resid_dge, design=edgerres_model)
resid_dge = estimateGLMTagwiseDisp(y=resid_dge, design=edgerres_model)
# fit GLM
resid_glm <- glmQLFit(y=resid_dge,design = edgerres_model)

edger_midp = GetMidPQuantileResiduals(resid_glm)


# ===== PCA =====
subtissue_cols = RColorBrewer::brewer.pal(n=length(unique(lung_annot$subtissue)),name="Dark2")
names(subtissue_cols) = unique(lung_annot$subtissue)
tmpCol = subtissue_cols["Inflamed Alveoli"]
subtissue_cols["Inflamed Alveoli"] = subtissue_cols["Normal Alveoli"]
subtissue_cols["Normal Alveoli"] = tmpCol


gVar = rowVars(edger_midp)
gVar = sort(gVar, decreasing=T)

pcout = prcomp(t(edger_midp[names(gVar)[1:1000],]),scale. = T, center = T)

ppca = ggplot(
  pcout$x %>% 
    as.data.frame() %>%
    rownames_to_column(var="Sample_ID") %>%
    left_join(
      y=lung_annot,
      by="Sample_ID"
    ),
  aes(x=PC1,y=PC2, color = subtissue, shape = subtissue)
) + geom_point() 


varExp = pcout$sdev**2
varExp = (varExp/sum(varExp))*100

ppcaFig = ggplot(
  pcout$x %>% 
    as.data.frame() %>%
    rownames_to_column(var="Sample_ID") %>%
    left_join(
      y=lung_annot,
      by="Sample_ID"
    ),
  aes(x=PC1,y=PC2, color = subtissue, shape = subtissue)
) + 
  geom_point(size=3) +
  scale_color_manual(values=subtissue_cols) +
  xlab(paste0("PC1 (",round(varExp[1],2),"%)")) +
  ylab(paste0("PC2 (",round(varExp[2],2),"%)")) +
  theme_cowplot() +
  theme(legend.title = element_blank())


ppcaFig


# ==== Inflamed vs Normal Alveoli ====
subsetAnnot = lung_annot %>%
  dplyr::filter(subtissue %in% c("Inflamed Alveoli", "Normal Alveoli"),Patient == 2)
subsetAnnot$subtissue = factor(subsetAnnot$subtissue )
subsetAnnot$subtissue = relevel(subsetAnnot$subtissue,ref="Normal Alveoli")

canonical_pathways = getGmt(
  "/data/work/Projects/BrScRNAseq/Plots/BroadPaper/CODE/data/c2.cp.v7.1.symbols.gmt"
)

cpLst2 = lapply(
  geneIds(canonical_pathways),
  function(x){
    x[x %in% rownames(lungDSP_cpm)]
  }
)


gsvaAlv <- GSVA::gsva(
  expr = edger_midp[,subsetAnnot$Sample_ID], 
  gset.idx.list = canonical_pathways, 
  method="ssgsea",
  kcdf="Gaussian"
)


# ===== Bubble Plot: Interferon-Gamma Signaling ======
# plot ssGSEA scores on the ROI positions 
slideDF = merge(
  x=data.frame(
    Sample_ID = subsetAnnot$Sample_ID,
    ssGSEA = gsvaAlv["REACTOME_INTERFERON_GAMMA_SIGNALING",subsetAnnot$Sample_ID]
  ),
  y=subsetAnnot,
  by="Sample_ID",
  all.x=TRUE
)

ggplot(
  data = slideDF,
  aes(x=x_coords,y=y_coords,size=ssGSEA,color=subtissue)
) + 
  geom_point(alpha=0.65) +
  scale_color_manual(values=subtissue_cols) +
  theme_cowplot() +
  labs(size = "ssGSEA Score", color = "") +
  #ggtitle("Interferon Gamma Signaling") +
  coord_fixed() +
  scale_x_reverse() +
  scale_y_reverse() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )


# ===== Heatmap: Interferon-Gamma Signaling ======
tmpAnnot = subsetAnnot %>% 
  dplyr::select(viral_score,subtissue) %>%
  dplyr::rename(`Viral Score` = viral_score,AOI = subtissue)
rownames(tmpAnnot) = subsetAnnot$Sample_ID
tmpAnnot$AOI = as.character(tmpAnnot$AOI)

tmpAnnot = tmpAnnot[order(tmpAnnot$AOI),]


pheatmap(
  mat = lungDSP_cpm[cpLst2[["REACTOME_INTERFERON_GAMMA_SIGNALING"]],rownames(tmpAnnot)],
  scale="row",
  cluster_cols = F,
  show_colnames = F,
  annotation_col = tmpAnnot %>% dplyr::select(AOI),
  annotation_colors = list(AOI = subtissue_cols[unique(tmpAnnot$AOI)]),
  fontsize_row = 5.75,
  fontsize = 9
)
