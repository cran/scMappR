#' Generate signature matrix 
#'
#' Convert a list of cell-type markers from FindMarkers in Seurat to a signature matrix defined by odds ratio and rank.
#'
#' Take a list of compiled differentially expressed genes from different cell-types, identify what the cell-types are using the Fisher's exact test, and then convert into a signature matrix for both the adjusted p-value and odds ratio.
#' 
#' @rdname generes_to_heatmap
#' @name generes_to_heatmap
#'
#' @param generes A list of cell-type markers with fold-changes and p-values (FindMarkers output in Seurat).
#' @param species The species of gene symbols, if not internal, "human" or "mouse".
#' @param naming_preference Likely cell-types given tissues (to be passed into human_mouse_ct_marker_enrich).
#' @param make_names Identify names of cell-type markers using the Fisher's exact test method (T/F).
#' @param internal If this function is pre-processing from Panglao (T/F).
#' @param rda_path Path to output direcotry, if toSave is true.
#'
#' @return List with the following elements:
#' \item{pVal}{A dataframe containing the signature matrix of ranks (-log10(Padj) * sign(fold-change)).}
#' \item{OR}{A dataframe containing the signature matrix of odds ratios.}
#' \item{cellname}{A vector of the cell-labels returned from the GSVA method.}
#' \item{topGenes}{the top 30 mos expressed genes in each cell-type.}
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
#' data(POA_example)
#'  POA_generes <- POA_example$POA_generes
#' signature <- generes_to_heatmap(POA_generes,species = -9, make_names = FALSE)
#' 
#' 
#' 
#' @export
#' 
generes_to_heatmap <- function(generes = generes,  
                               species = "human",  
                               naming_preference = -9, 
                               rda_path = "", 
                               make_names = TRUE, 
                               internal = FALSE) {
  # take an list of compiled DEGs from different cell types, identify what the cell-types are using the fisher's exact method, and then convert into a signature matrix for both the adjusted p-value and odds ratio
  # Args:
  # generes: A list of cell-type markers where each element of the list has the gene symbol, p_adj, and fold-change: the output of the FindMarkers function in Seurat
  # colnames are asfollows
  #  p_val_adj avg_logFC
  # species: human, mouse, or -9 if internal
  # naming_preference: if you have an idea of what cell-types should be in the dataset, then using one of the arguments in the 
  # get_naming_preference() funciton will increase the rank that cell-types falling into that category
  # make_names: If you ant the cell-types to be named or not
  # internal: If this is part of our internal pipeline, we apend the mouse and human symbol onto the final matrix so that we can easily tell afterwards
  # Returns:
  # A list containing a signature matrix by rank := -1*log10(Pfdr) and by fold-change (only increasing). 
  # additionally it returns the top (up to) 30 CT markers for each cell-type, as well as the name of each cell-type (from the signature methods method)
  
  if(!(is.list(generes))) {
    stop("generes object must be a list.")
  }
  
  if(!(species %in% c("human", "mouse"))) {
    if(species != -9) {
      stop("species is not 'human' 'mouse' or '-9' (case sensitive), please try again with this filled.")
    }
  }
  naming_preferences <- c("brain", "epithelial", "endothelial", "blood", "connective","eye", "epidermis", "Digestive", "Immune", "pancreas", "liver", "reproductive", "kidney", "respiratory") 
  if(!naming_preference %in% naming_preferences) {
    if(naming_preference != -9) {
      message("Naming preference options")
      message(naming_preferences)
    stop("Naming preferences not in options (case sensitive) and isn't a non-choice (-9), please try again.")
    }
  }
  
  if(!(is.character(rda_path))) {
    stop("rda_path must be of class character")
  }
  if(all(is.logical(make_names), is.logical(internal))[1] == FALSE) {
    stop("make_names and internal must be class logical (TRUE/FALSE)")
  }
  
  if(species == -9) {
    # if it's internal and symbols and ENSBML are attached
    for(i in 1:length(generes)) {
      generes[[i]] <- generes[[i]][stats::complete.cases(generes[[i]]),]
      names1 <- get_gene_symbol(generes[[i]])
      rownames(generes[[i]]) <- names1$rowname
    }
    species <- names1$species
  }
  topGenes <- topgenes_extract(generes) # take the top 30 genes
  if(make_names == TRUE) {
    cell_name <- human_mouse_ct_marker_enrich(topGenes, theSpecies = species, cell_marker_path = rda_path, naming_preference = naming_preference) 
    # get the names of each of the cell types
    cell_name <- cell_name$cellTypes # attach the appropriate cell names
  } else {
    cell_name = colnames(generes)
  }
  genes <- c()
  for(i in 1:length(generes)) {
    genes <- c(genes,rownames(generes[[i]]))
  }
  # generate matrix with ranks or odds ratio
  # only genes that are preferentially expressed in at min one cell type
  # is removed
  genes_uni <- unique(genes)
  if(internal == TRUE) {
    if(species == "human") {
      genes_uni <- paste0(genes_uni, "-ENSG00000000000.0")
      for(i in 1:length(generes)) {
        rownames(generes[[i]]) <- paste0(rownames(generes[[i]]),"-ENSG00000000000.0")
      }
      
    }
    if(species == "mouse") {
      genes_uni <- paste0(genes_uni, "-ENSMUSG00000000000.0")
      for(i in 1:length(generes)) {
        rownames(generes[[i]]) <- paste0(rownames(generes[[i]]),"-ENSMUSG00000000000.0")
      }
    }
  }
  
  if(is.null(generes[[1]]$avg_logFC)) {
    message("Seurat V4 or later was used to identify cell-type markers, adding 'avg_logFC' column. It has the same data as avg_log2FC but is compatible with downstream functions in scMappR.")
    for(z in 1:length(generes))
      generes[[z]]$avg_logFC <- generes[[z]]$avg_log2FC
  }
  
  # generating the signature matrix developed by the rank of p-vlaues
  scmappr <- matrix(0, nrow = length(genes_uni), ncol = length(names(generes))) 
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- -1*log10(generes[[i]]$p_val_adj) * sign(generes[[i]]$avg_logFC)
    # compute the padj and fold change into the rank
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_rank_mat_t <- scmappr
  wilcoxon_rank_mat_t[is.infinite(wilcoxon_rank_mat_t) & wilcoxon_rank_mat_t < 0] <- min(wilcoxon_rank_mat_t[is.finite(wilcoxon_rank_mat_t)])
  wilcoxon_rank_mat_t[is.infinite(wilcoxon_rank_mat_t) & wilcoxon_rank_mat_t > 0] <- max(wilcoxon_rank_mat_t[is.finite(wilcoxon_rank_mat_t)])
  
  
  # generating the signature matrix developed by the rank of fold changes
  scmappr <- matrix(-300, nrow = length(genes_uni), ncol = length(names(generes)))
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- generes[[i]]$avg_logFC
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_rank_mat_or <- scmappr
  wilcoxon_rank_mat_or[is.infinite(wilcoxon_rank_mat_or) & wilcoxon_rank_mat_or < 0] <- min(wilcoxon_rank_mat_or[is.finite(wilcoxon_rank_mat_or)])
  wilcoxon_rank_mat_or[is.infinite(wilcoxon_rank_mat_or) & wilcoxon_rank_mat_or > 0] <- max(wilcoxon_rank_mat_or[is.finite(wilcoxon_rank_mat_or)])
  wilcoxon_rank_mat_or <- 2^wilcoxon_rank_mat_or # instead of the log2Fc as the output, ouput the odds ratio (only positive) that the gene is in that cell-type
  l <- list(pVal = wilcoxon_rank_mat_t, OR = wilcoxon_rank_mat_or, cellname = cell_name, topGenes = topGenes)
  # return the pVal signature matrix, odds ratio signautre matrix, cell names for that list, and topGenes
  return(l)
}

