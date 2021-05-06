### Convert scanpy to formats for R

#conda activate scanpy
#import anndata as ad
#adata=ad.read_h5ad("~/broad_kidney_review.h5ad")
# import scipy.io
# scipy.io.mmwrite("~/covid_kidney_combined.mtx",adata.layers['winsorized'])
#  adata.var.to_csv("~/covid_kidney_combined_var.csv")
# adata.obs.to_csv("~/covid_kidney_combined_obs.csv")


library(Seurat)
library(readr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_cowplot())
library(harmony)
library(magrittr)
library(Hmisc)
library(hdf5r)


obj=readMM("~/covid_kidney_combined.mtx")
obj=t(obj)
obs=read_csv("~/covid_kidney_combined_obs.csv")
var=read_csv("~/covid_kidney_combined_var.csv")
rownames(obj)=var$featurekey
colnames(obj)=obs$barcodes

seurat=CreateSeuratObject(obj,project="covid_kidney")
df=data.frame(sample=obs$sample)
rownames(df)=rownames(seurat@meta.data)
seurat=AddMetaData(seurat,metadata=df)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

seurat=NormalizeData(seurat,verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = 20, verbose = FALSE)

seurat <- seurat %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

seurat_markers=FindAllMarkers(seurat,min.pct=0.25)

##################subclustering#############################

###### immune/CD45+ ######

kidney_imm = subset(seurat, ident=c(13,14,20))
kidney_imm <- kidney_imm %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = kidney_imm@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

kidney_imm_markers=FindAllMarkers(kidney_imm)

###### EC/PECAM1+ ######

kidney_ec = subset(seurat, ident=c(6,15))
kidney_ec <- kidney_ec %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = kidney_ec@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

kidney_ec_markers=FindAllMarkers(kidney_ec)

###### PCT ######

kidney_PCT = subset(seurat, ident=c(0,1,2,21,22))
kidney_PCT <- kidney_PCT %>% FindVariableFeatures(selPCTtion.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = kidney_PCT@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

kidney_PCT_markers=FindAllMarkers(kidney_PCT)

###### mes ######

kidney_mes = subset(seurat, ident=c(8,10))
kidney_mes <- kidney_mes %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = kidney_mes@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

kidney_mes_markers=FindAllMarkers(kidney_mes)

### Annotations ######
file=read_tsv(paste0("kidney_celltypes.txt"))
seurat <- AddMetaData(seurat,metadata=file)
Idents(seurat) <- "annotation"


##########FIGURES#############


df=read_tsv(paste0("ED10def_kidney_sourcedata.tsv"))

######ED10D#########
df%>%ggplot(aes(x=UMAP_1,y=UMAP_2,col=celltype))+geom_point(size=1)+scale_colour_manual(values=c("#a20a0a","#ffcac2","#f6eec9","#a4b787","#cbaf87","#ff9a76","#9ad3bc","#79d70f","#bbbfca","#1f6f8b","#790c5a","#e5c5b5","#d2f5e3","#9ab3f5","#213e3b","#f0a500",c25[-6]))+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))

######ED10E#########
df1%>%ggplot(aes(x=sample,y=1,fill=celltype))+geom_col(position="fill")+scale_fill_manual(values=c("#a20a0a","#ffcac2","#f6eec9","#a4b787","#cbaf87","#ff9a76","#9ad3bc","#79d70f","#bbbfca","#1f6f8b","#790c5a","#e5c5b5","#d2f5e3","#9ab3f5","#213e3b","#f0a500",c25[-6]))+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))


######ED10F#########
donors=read_tsv(paste0("kidney_donor_colors.txt"))
df%>%ggplot(aes(x=UMAP1,y=UMAP2,col=Donor))+geom_point(size=1)+scale_colour_manual(values=donors$color)+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))

######ED11A#########
DotPlot(seurat,features=c("ACE2","TMPRSS2","CTSL"))+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))







