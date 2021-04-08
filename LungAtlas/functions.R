# A function that subsamples the single-cell dataset while maintaining
# representation from multiple cell types and performs proportion estimation 
# using MuSiC
#
# input.bulk.eset bulk ExpressionSet Object
# input.sc.eset: single-cell ExpressionSet Object
# n.cells: number of cells to downsample to (default: 10k)
subsample_and_estimate_proportion <- function(input.bulk.eset, input.sc.eset, n.cells=10000) {
  # set samplename (MuSiC requirement)
  pData(input.sc.eset)$samplename <- input.sc.eset$barcodes
  
  # Perform stratified sampling
  tmp <- as.data.frame(pData(input.sc.eset))
  cell.sample.idx <- which(pData(input.sc.eset)$barcodes %in% splitstackshape::stratified(tmp, group='Cluster',size=n.cells/nlevels(tmp$Cluster) )$barcodes)
  sc.eset.sample <- input.sc.eset[,cell.sample.idx]
  
  # Perform deconvolution
  est.prop.all <- music_prop(bulk.eset=input.bulk.eset,
                             sc.eset=sc.eset.sample,
                             clusters='Cluster',
                             samples='samplename')
  
  # Return results
  est.prop.all
}
