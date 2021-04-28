# covid19-autopsy

This repository includes the code to process and analyze the data from 17 lung autopsy donors that succumed to COVID-19. This study was titled "COVID-19 tissue atlases reveal SARS-COV-2 pathology and cellular targets" and published in Nature on April 2021.

The code is organized in directories by the analysis that were done:

**DataPreprocessing** - Pre-processing done on lung sc-RNAseq and sn-RNAseq of lung samples, including Cell Bender for ambient RNA removal, doublet removal and quality control.

**ExtendedFigure9** - Inflammation-specific programs in spatial lung samples. Scripts to generate results in Extended Figure 9.

**LungAnnotationPrediction** - Code to predict cell type annotations in the lung.

**LungBulk** - Analysis of bulk RNAseq data of lung samples including deconvolution.

**LungDE** - Differential expression analysis of COVID lung samples to healthy lung datasets from prior published studies.

**OtherTissues** - Initial analysis of other tissues that were collected in this study.

**SpatialAnalysis** - Spatial analysis of lung samples. This directory includes a detailed README file to explain its content.

**Sub-Clustering** - Manual annotaitons of lung cell types and sub-clusterings of the following main lineages: endotheila cells, B+Plasma cells, T+NK cells and myeloid cells.

**ViralAnalysis** - Viral enrichment analysis for lung and heart, and the identification of cell types that are enriched in viral RNA.
