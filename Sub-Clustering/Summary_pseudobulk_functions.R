require(stringr)
require(limma)
require(Matrix)
require(dplyr)
require(magrittr)
require(data.table)
require(glue)
require(pbapply)
require(reticulate)
require(SingleCellExperiment)
require(Seurat)
require(openxlsx)
require(ggplot2)


# A is a dgCMatrix of Genes(rows) X leiden(columns)
do_log2cpm <- function(A, total = 1e4) 
{
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- total * A@x
  A@x <- log2(A@x + 1)
  return(A)
}


CalcPseudoBulk = function(data, genes, meta, title, with_patient_cor)
{
  # Organizing data
  stopifnot(nrow(genes)==nrow(data))
  stopifnot(nrow(meta)==ncol(data))
  
  rownames(data) = genes[,1]
  colnames(data) = rownames(meta)
  
  if (with_patient_cor=="True")
    y <- with(meta, model.matrix(~ 0 + factor(leiden):factor(patient)))
  else
    y <- with(meta, model.matrix(~ 0 + factor(leiden)))
  y <- as(y, "dgCMatrix")
  pb <- as(data %*% y, "dgCMatrix") #Summarize all cells per cluster
  pb <- do_log2cpm(pb, median(colSums(pb))) # normalize

  if (with_patient_cor=="True")
  {
      pb_meta <- str_split_fixed(colnames(pb), ":", 2)
      colnames(pb_meta) <- c("leiden", "patient")
      pb_meta <- as_tibble(pb_meta)
      pb_meta %<>%
        mutate(
          leiden = str_replace(leiden, "factor\\(leiden\\)", ""),
          patient = str_replace(patient, "factor\\(patient\\)", "")
        )
   }
   else
       pb_meta=meta
   stopifnot(nrow(pb_meta) == ncol(pb))

  ################################
  # All versus All (AVA)
  # Test all pairs of clusters
  ################################
  print("****** All vs all analysis *******")
  pb_meta$x <- factor(pb_meta$leiden)
  des1 <- with(pb_meta, model.matrix(~ 0 + x)) # Design matrix split by leiden cluster+patient
  fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5,]), design = des1) # Testing only rows with mean expression above a threshold (to help multiple hypothesis testing)
  fit1 <- eBayes(fit1)
  fit1$genes <- rownames(fit1$coefficients)
  
  ## (a) Comparing Pairwise clusters
  cluster_pairs <- t(combn(levels(pb_meta$x), 2))
  cont <- makeContrasts(contrasts = lapply(seq(nrow(cluster_pairs)), function(i) 
  {
    glue("x{cluster_pairs[i,1]} - x{cluster_pairs[i,2]}")
  }), levels = des1)
  
  colnames(cont) <- str_replace(colnames(cont), " - ", "vs")
  fit2 <- contrasts.fit(fit1, cont)
  fit2 <- eBayes(fit2)
  
  de_ava <- rbindlist(lapply(colnames(cont), function(this_coef) 
  {
    x <- topTable(fit2, coef = this_coef, number = nrow(fit1$coefficients))
    this_coef <- str_replace_all(this_coef, "x", "")
    this_coef <- str_replace(this_coef, "vs", " vs ")
    x$coef <- this_coef
    x
  }))
  
  de_ava_file <- paste0(title,"_de_AVA.xlsx")
  print(glue("Writing {de_ava_file}"))
  all_comps = unname(unlist(as.list(unique(de_ava[,"coef"]))))
  tabs_list=list()
  for(comp in all_comps)
    tabs_list[[comp]] = subset(de_ava, coef==comp)
  write.xlsx(tabs_list, file = de_ava_file)

      # Specific to infected vs. non-infected PB
      pdf(paste0("VolcanoPlot_",title,".pdf"))
      DE_df = subset(de_ava,coef=="False vs True")
      DE_df[,"min_log_adj_pval"] = -log(DE_df[,"adj.P.Val"])
      print(ggplot(DE_df, aes(logFC, min_log_adj_pval)) + geom_point() + geom_vline(xintercept = c(-2,2), col="red") +
         geom_hline(yintercept = c(2), col="red") + ggtitle(paste0("DE genes - ", title)))
      dev.off()


  ################################
  # One versus All (OVA)
  ################################
  print("****** One vs all analysis *******")
  
  de_ova <- rbindlist(pblapply(sort(unique(pb_meta$leiden)), function(this_cluster) 
  {
    pb_meta$x <- pb_meta$leiden == this_cluster
    des1 <- with(pb_meta, model.matrix(~ x))
    fit1 <- lmFit(object = as.matrix(pb[rowMeans(pb) > 0.5,]), design = des1)
    fit1 <- eBayes(fit1)
    fit1$genes <- rownames(fit1$coefficients)
    res <- topTable(fit1, coef = 2, number = 1e6)
    res$coef <- sprintf("%s vs all", this_cluster)
    return(res)
  }))
  de_ova_file <- glue(title,"_de_OVA.xlsx")
  print(glue("Writing {de_ova_file}"))
  all_comps = unname(unlist(as.list(unique(de_ova[,"coef"]))))
  tabs_list=list()
  for(comp in all_comps)
    tabs_list[[comp]] = subset(de_ova, coef==comp)
  write.xlsx(tabs_list, file = de_ova_file)

}
