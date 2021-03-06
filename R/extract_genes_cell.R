#' Extract Markers
#'
#' Extracting cell-type markers from a signature matrix.
#'
#' This function takes a signature matrix and 
#' extracts cell-type markers above a p-value or fold-change threshold.
#'
#'
#' @rdname extract_genes_cell
#' @name extract_genes_cell
#'
#' @param geneHeat The heatmap of ranks from your scRNA-seq dataset with your genes subsetted.
#' @param cellTypes The cell-types that you're interested in extracting. They need to be colnames (not case sensitive).
#' @param val How associated a gene is with a particular cell type to include in your list - default is slightly associated.
#' @param isMax If you are taking the single best CT marker (T/F) -- TRUE not recommended.
#' @param isPvalue If the signature matrix is raw p-value (T/F) -- TRUE not recommended.
#' 
#'
#' @return \code{extract_genes_cell} A vector of genes above the threshold for each sample. \cr
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text geom_bar theme_classic xlab ylab scale_fill_manual element_line
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test dist shapiro.test mad
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#' @importFrom reshape melt
#'
#' @examples
#' 
#' 
#' data(POA_example)
#' Signature <- POA_example$POA_Rank_signature
#' RowName <- get_gene_symbol(Signature)
#' rownames(Signature) <-RowName$rowname
#' # extract genes with a -log10(Padj > 1)
#' Signat <- extract_genes_cell(Signature)
#' 
#' 
#' 
#' @export
#' 
extract_genes_cell <- function(geneHeat, cellTypes = "ALL", val = 1, isMax = FALSE, isPvalue = FALSE) {

# This function takes a signature matrix and extracts cell-type markers above a p-value or fold-change threshold. 

  #Args:   
    # geneHeat: the heatmap of ranks from your scRNA-seq dataset with your genes subsetted
    # cellTypes: The cell-types that you're interested in extracting. They need to be colnames (not-case sensitive)
    # val: how associated a gene is with a particualr cell type to include in your list - default is slightly above 1
    # isMax true or false, if you want to sort genes into the CT marker that represents theme best
    # isPvalue: if the signature matrix is simply a raw P-value: not reccomended
  #Returns: a list of genes above the threshold for each sample.

  geneHeat_class <- class(geneHeat)[1] %in% c("data.frame", "matrix")
  if(geneHeat_class[1] == FALSE ) {
    stop("geneHeat should be a data.frame or matrix of ranks from your scRNA-seq dataset.")
  }
  
  if(!is.character(cellTypes)) {
    stop("cellTypes should be of class character -- either column names of cell-types to include or 'ALL' -- ALL is reccomended.")
  }
  
  if(!is.numeric(val)) {
    stop("val must be of class numeric.")
  }
  if(all(is.logical(isMax), is.logical(isPvalue))[1] == FALSE) {
    stop("isMax and isPvalue must be of class logical (TRUE/FALSE).")
  }
  
colnames(geneHeat) <- toupper(colnames(geneHeat))
cellTypes <- toupper(cellTypes)

geneHeat <- geneHeat[rowSums(geneHeat) > 1,] # extract genes with any CT specificity

if(any(is.character(geneHeat), is.numeric(geneHeat))[1] ) {
  if((!is.matrix(geneHeat))[1]) {
  geneHeat <- data.frame(t(geneHeat))
  }
}  

genes_extracted <- list()
if(cellTypes == "ALL") { # if you extract all cell-type just make it your output
  cellTypes <- colnames(geneHeat)
}

if(isMax == TRUE) {  #take the single top gene
  # this way, you take the cell type that marks each gene the best.
  # no gene is double counted into multiple cell types.
  warning("Extracting the single most expressed marker in each cell-type. Not reccomeneded.")
  cn <- colnames(geneHeat)[max.col(geneHeat,"first")] # the cell type with the most signficant
  names(cn) <- rownames(geneHeat)
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- names(cn)[which(cn == colnames(geneHeat)[i])]
    names(genes_extracted)[i] <- colnames(geneHeat)[i]
  }
  return(genes_extracted) 
} else { # take top genes past a threshold
  # a gene can be counted into multiple cell types and you're just saying that a 
  # gene is in it if it's p-value passes a certain threshold -- This is probably
  # what fits the premise of scMappR slightly better
  if(isPvalue == TRUE) {
    val = -1*log10(val)
  } else {
    val = val
  }
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- rownames(geneHeat)[geneHeat[,i] > val]
    names(genes_extracted)[i] <- names(genes_extracted)[i] <- colnames(geneHeat)[i]
    
  }
  return(genes_extracted)
}
}
