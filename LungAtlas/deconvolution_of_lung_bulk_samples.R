# Deconvolution of Bulk Lung Samples

# Load Libraries
library(data.table)
library(Matrix)
library(MuSiC)
library(xbioc)
library(tidyverse)
library(tidyverse)
library(cowplot)
library(fastSave)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(nbHelpers)
library(ggpattern)

# Load auxilary codes
source('functions.R')

# Samples that were present in the input data but are not part of bulk lung samples
exclude.samples <- c("02_P334354_S106_R01","02_P166169_S027_R01","04_P054921_S029_R01")
exclude.samples.donor <- c("D6","D3")

# Location of input datasets
bulk.input.path <- "path_to_bulk_input_data"
sc.input.path <- "path_to_single_cell_input_data"

# Load bulk datasets
bulk.counts <- fread(paste0(bulk.input.path, 'rsem.genes.counts.matrix'))
bulk.matrix <- as.matrix(bulk.counts[,-1])
rownames(bulk.matrix) <- bulk.counts$V1

# Create ExpressionSet object
bulk.eset <- ExpressionSet(
  assayData = as.matrix(bulk.matrix),
  phenoData = AnnotatedDataFrame(data.frame(row.names=colnames(bulk.matrix),samplename=colnames(bulk.matrix))),
  featureData = AnnotatedDataFrame(data.frame(row.names=rownames(bulk.matrix),samplenames=rownames(bulk.matrix)))
)

# Inspect ExpressionSet object
head(pData(bulk.eset))
head(fData(bulk.eset))
colnames(bulk.eset)
dim(bulk.eset)

# Remove samples not used here from the bulk
bulk.eset <- bulk.eset[,!(colnames(bulk.eset) %in% exclude.samples)]

# Load Single Cell dataset
sc.raw <- readMM(paste0(sc.input.path,'/raw.mtx'))
cell_metadata <- read.csv(paste0(sc.input.path,'cell_metadata.csv'))
gene_metadata <- read.csv(paste0(sc.input.path,'gene_metadata.csv'))
rownames(sc.raw) <- cell_metadata$barcodes
colnames(sc.raw) <- gene_metadata$featurekey
head(cell_metadata)
rownames(cell_metadata) <- cell_metadata$barcodes

# Make Single Cell ExpressionSet
sc.eset <- ExpressionSet(as.matrix(t(sc.raw)), 
                         phenoData = AnnotatedDataFrame(cell_metadata))


#############################################
## Prepare the data for deconvolution
#############################################

# Remove cell types that represent unclassified cells
discard.cell.types <- c('(sample specific)','Doublets')
sc.eset.clean <- sc.eset[,!as.character(pData(sc.eset)$Cluster) %in% discard.cell.types]

# Update Cluster variable levels
pData(sc.eset.clean)$Cluster <- as.factor(as.character(pData(sc.eset.clean)$Cluster))


# Perform 10 independent subsamplings of the data followed by proportion estimation
prop.estimate.perm <- lapply(1:10, function(i) { 
  # perform garbage collection
  gc()
  # run iteration
  subsample_and_estimate_proportion(input.bulk.eset=bulk.eset, 
                                    input.sc.eset=sc.eset.clean,
                                    n.cells=10000)
})

# Merge the output data into a single table and append iteration nuber
perm_melt_data <- do.call(rbind, lapply(seq_along(prop.estimate.perm), function(i){
  prop_melt <- reshape2::melt(prop.estimate.perm[[i]]$Est.prop.weighted)
  prop_melt$iter <- i
  prop_melt
}))

# Fix the header names
head(perm_melt_data)
colnames(perm_melt_data) <- c("sample","celltype","pc","iter")

# Get an ordering for the celltypes based on mean abundance
perm_melt_data %>% 
  filter(!(sample %in% exclude.samples)) %>% 
  group_by(celltype) %>% 
  summarise(avg=mean(pc)) %>%
  arrange(desc(avg)) %>% pull(celltype) -> ct.order

ct.order <- as.character(ct.order)
ct.order

# Generate a color palette
pal <- brewer.pal(name='Set3',n=length(ct.order))
names(pal) <- ct.order

#################################
## Fix identifiers for final plot
#################################

# Update identifiers
perm_melt_data$celltype <- factor(perm_melt_data$celltype, levels = ct.order)

# Create short sample identifiers for matching with id_map
perm_melt_data$sample_short <- unlist(lapply(strsplit((as.character(perm_melt_data$sample)), '_', '['), function(s) {paste(s[1:2], collapse='-')}))

## Load Id map
id_map <- read.table('id_map.txt',header=T,stringsAsFactors = FALSE)

# Set Donor_ID variable
perm_melt_data$Donor_ID <- id_map$Donor_ID[match( perm_melt_data$sample_short,id_map$Broad_ID )]
perm_melt_data$Donor_ID <- factor(perm_melt_data$Donor_ID)
perm_melt_data$Donor_ID <- factor(perm_melt_data$Donor_ID,levels=paste0('D',sort(as.numeric(gsub('D','',levels(perm_melt_data$Donor_ID))))))

#######################################################
# Generate estimates with increasing number of cells
#######################################################

# The purpose of this section is to demonstrate that there is no 
# systematic bias introduced when the number of cells used is small

# Generate cell counts
prop.estimate.perm.increasing <- lapply(rep(seq(5000,10000,500),each=2), function(n.cells) { 
  subsample_and_estimate_proportion(input.bulk.eset=bulk.eset, input.sc.eset=sc.eset.clean, n.cells=n.cells)
})
names(prop.estimate.perm.increasing) <- rep(seq(5000,10000,500),each=2)

# Merge the data into a single data.frame
z <- lapply(prop.estimate.perm.increasing, function(t) {
  reshape2::melt(t$Est.prop.weighted)
})
zz <- lapply(seq_along(z), function(i) {
  entry <- z[[i]]
  ncells <- as.numeric(names(z)[i])
  entry$ncells <- c(ncells)
  entry
})
zzz <- do.call(rbind, zz)

# Check the output
head(zzz)

# Fix columns
zzz$sample_short <- unlist(lapply(strsplit((as.character(zzz$Var1)), '_', '['), function(s) {paste(s[1:2], collapse='-')}))
zzz$Donor_ID <- id_map$Donor_ID[match( zzz$sample_short,id_map$Broad_ID )]
zzz$Donor_ID <- factor(zzz$Donor_ID)
zzz$Donor_ID <- factor(zzz$Donor_ID,levels=paste0('D',sort(as.numeric(gsub('D','',levels(zzz$Donor_ID))))))

# Generate Plot "robust_to_subsampling"
zzz %>% 
  filter(!(Var1 %in% exclude.samples)) %>%
  ggplot( aes(x=ncells,y=value,color=Var2)) + 
  geom_point() + facet_wrap(~Donor_ID) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text()) +
  xlab("Number of cells in reference") +
  ylab("% of total") + labs(color="Cell types") +
  scale_color_manual(values=pal)

# Save as pdf
ggsave("robust_to_subsampling.pdf",width=10,height=10)

# Get the single-cell metadata
sc.metadata <- as.data.frame(pData(sc.eset.clean))
sc.metadata[!sc.metadata$Cluster %in% discard.cell.types,] -> df1

# Get cluster/sample counts
df1 %>% group_by(sample,Cluster) %>% count() -> df2
dcast(df2, Cluster ~ sample) -> df3
rownames(df3) <- df3$Cluster
df3 <- df3[,-1]

# Remove NAs
df3[is.na(df3)] <- 0

# Calculate proportions
sweep(df3, 2,STATS = apply(df3, 2, sum), FUN = '/') -> proportions

# convert proportions to long format
prop.m <- melt(as.matrix(proportions))

# Fix column names
colnames(prop.m) <- c('Celltype','Sample','pc')

# Update Celltype order
prop.m$Celltype <- factor(prop.m$Celltype,levels = ct.order) 

# Amend sample identifiers
prop.m$Sample_short <- gsub('-S.*','',prop.m$Sample)
prop.m$Donor.id <- id_map$Donor_ID[match(prop.m$Sample_short, id_map$Broad_ID)]

# Generate estimated proportions data frame
perm_melt_data %>% 
  filter(!(sample %in% exclude.samples)) -> estimated

# Generate single-cell proportions data frame
prop.m %>% 
  filter(!(as.character(Sample) %in% gsub('_','-',exclude.samples))) %>%
  group_by(Celltype,Donor.id) %>% 
  summarize(pc=mean(pc)) %>%  
  mutate(Donor.id = factor(Donor.id, levels=(paste0('D',c(1:17))))) -> sc.prop

# Make adjustments to bulk proportions data set format
estimated %>% 
  select(celltype, Donor_ID, pc) %>% 
  mutate(type=c('estimated')) -> merge1

# Make adjustment to single-cell proportions data set format
sc.prop %>% 
  mutate(type='actual') -> merge2

# Convert both to data.frame
merge1 <- as.data.frame(merge1)
merge2 <- as.data.frame(merge2)

# Make column names identical
colnames(merge2) <- colnames(merge1)

# Merge the two datasets (estimated and single-cell proportions)
merged <- rbind(merge1, merge2)

# Inspect
head(merged)

# Generate a unified palette
tmp_pal_1 <- pal
tmp_pal_2 <- pal
names(tmp_pal_1) <- paste0(names(tmp_pal_1),'_actual')
names(tmp_pal_2) <- paste0(names(tmp_pal_2),'_estimated')
pal_new <- c(tmp_pal_1,tmp_pal_2)

# Generate mean and sd of the percent values 
merged %>% 
  group_by(celltype, Donor_ID, type) %>%
  summarize(pc0 = mean(pc), pc_sd = sd(pc) ) %>% 
  mutate(celltype_type = paste(celltype,type,sep='_')) -> tmp1

# inspect
head(tmp1)
  
# order the celltype_type factor
tmp1 %>% ungroup() %>% 
  select(celltype, type, celltype_type) %>% 
  unique() %>% 
  arrange(celltype, type) %>%
  pull(celltype_type) -> factor_order

# update columns
tmp1 %>% 
  mutate(celltype_type = factor(celltype_type, levels=factor_order)) -> tmp2

# Set any NAs to 0
tmp2$pc_sd[is.na(tmp2$pc_sd)] <- 0
    
# Generate the plot of celltype proportions for the single-cell and 
# estimated data
tmp2 %>%
  ggplot(aes(x=celltype_type,y=pc0,pattern=type,fill=celltype_type)) + 
  geom_bar_pattern(stat='identity', position='dodge', 
                   pattern_angle=45, pattern_density=0.1,
                   pattern_spacing=0.025,color='black') + 
  geom_errorbar(aes(ymin=pc0-pc_sd, ymax=pc0+pc_sd)) +
  facet_wrap(~Donor_ID) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Cell type') + 
  ylab("% composition") +
  scale_pattern_manual(values = c(estimated = "stripe", actual = "none"))  +
  scale_fill_manual(values=pal_new) +
  theme(legend.position = 'none')
  
# Save the final plot
ggsave('side_by_side_estimated_sc_errorbars.pdf',w=10,h=10)  
