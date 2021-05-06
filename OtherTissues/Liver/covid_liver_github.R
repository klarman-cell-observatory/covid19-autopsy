


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


####### LIVER #######
# import scipy.io
# import scanpy as sc
# adata=sc.read("~/broad_liver_prelim.h5ad")
# scipy.io.mmwrite("~/covid_liver_v.mtx",adata.layers['winsorized'])
# adata.var.to_csv("~/covid_liver_var.csv")
# adata.obs.to_csv("~/covid_liver_obs.csv")


obj=readMM("~/covid_liver_v.mtx")
obj=t(obj)
adata.var.to_csv("~/covid_liver_var.csv")
obs=read_csv("~/covid_liver_obs.csv")
var=read_csv("~/covid_liver_var.csv")
rownames(obj)=var$featurekey
colnames(obj)=obs$barcodes

liver=CreateSeuratObject(obj,project="covid_liver")

df=data.frame(sample=obs$sample)
rownames(df)=rownames(liver@meta.data)

liver=AddMetaData(liver,metadata=df)
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")
liver=NormalizeData(liver,verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = liver@var.genes, npcs = 20, verbose = FALSE)


liver <- liver %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

liver_markers=FindAllMarkers(liver,min.pct=0.25)


###### liver EC/PECAM1+ ######

liver_ec=subset(liver,ident=c(2,14,17,22,23))

liver_ec <- liver_ec %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat_ec@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

liver_ec_markers=FindAllMarkers(liver_ec)


###### liver epi ######


liver_epi=subset(liver,ident=c(0,1,3,4,5,8,11,13,20,21))

liver_epi <- liver_epi %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat_epi@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

liver_epi_markers=FindAllMarkers(liver_epi)


###### liver immune ######

liver_imm=subset(liver,ident=c(6,9,12,18,19))

liver_imm <- liver_imm %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat_imm@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

liver_imm_markers=FindAllMarkers(liver_imm)


###### mes ########

liver_mes=subset(liver,ident=c(7,16))

liver_mes <- liver_mes %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat_mes@var.genes, npcs = 20, verbose = FALSE) %>% RunHarmony("sample", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

liver_mes_markers=FindAllMarkers(liver_mes)


###### Liver figures ######

file=read_tsv(paste0(wd,"liver_celltypes.txt"))

liver <- AddMetaData(liver,metadata=df1_sub)

Idents(liver) <- "celltype"

########## Figures#########

df=read_tsv(paste0("ED10ghi_liver_sourcedata.tsv"))

##### ED10G#########
df%>%ggplot(aes(x=UMAP_1,y=UMAP_2,col=celltype))+geom_point(size=1)+scale_colour_manual(values=c("#a20a0a","#ffcac2","#f6eec9","#a4b787","#cbaf87","#ff9a76","#9ad3bc","#79d70f","#bbbfca","#1f6f8b","#790c5a","#e5c5b5","#d2f5e3","#9ab3f5","#213e3b","#f0a500",c25[-6]))+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))

##### ED10H#########
df%>%ggplot(aes(x=sample,y=1,fill=celltype))+geom_col(position="fill")+scale_fill_manual(values=c("#a20a0a","#ffcac2","#f6eec9","#a4b787","#cbaf87","#ff9a76","#9ad3bc","#79d70f","#bbbfca","#1f6f8b","#790c5a","#e5c5b5","#d2f5e3","#9ab3f5","#213e3b","#f0a500",c25[-6]))+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))

##### ED10I#########
donors=read_tsv(paste0("liver_donor_colors.txt"))
df%>%ggplot(aes(x=UMAP1,y=UMAP2,col=Donor))+geom_point(size=1)+scale_colour_manual(values=donors$color)+ guides(colour = guide_legend(ncol=2,override.aes = list(size=5),title=""))


######ED11B#########
DotPlot(liver,features=c("ACE2","TMPRSS2","CTSL"))+theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))








