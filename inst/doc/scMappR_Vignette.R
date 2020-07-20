## ----install_developter, eval=FALSE-------------------------------------------
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  
#  BiocManager::install("pcaMethods")
#  BiocManager::install("GSVA")
#  
#  devtools::install_github("DustinSokolowski/scMappR")
#  
#  
#  

## ----install_cran, eval=FALSE-------------------------------------------------
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  if (!requireNamespace("devtools", quietly = TRUE))
#      install.packages("devtools")
#  
#  BiocManager::install("pcaMethods")
#  BiocManager::install("GSVA")
#  
#  install.packages("scMappR")
#  

## ----library_scMappR, warning=FALSE, echo = FALSE-----------------------------

library(scMappR)


## ----scMappR_internal_example, eval = FALSE-----------------------------------
#  
#  data(POA_example) # region to preoptic area
#  
#  Signature <- POA_example$POA_Rank_signature # signature matrix
#  
#  rowname <- get_gene_symbol(Signature) # get signature
#  
#  rownames(Signature) <- rowname$rowname
#  
#  genes <- rownames(Signature)[1:60]
#  
#  rda_path1 = "" # data directory (if it exists)
#  
#  # Identify tissues available for tissue_scMappR_internal
#  data(scMappR_tissues)
#  
#  "Hypothalamus" %in% toupper(scMappR_tissues)
#  
#  internal <- tissue_scMappR_internal(genes, "mouse", output_directory = "scMappR_Test_Internal",
#  tissue = "hypothalamus", rda_path = rda_path1, toSave = TRUE, path = tempdir())
#  
#  

## ----scMappR_custom_example, eval = FALSE-------------------------------------
#  
#  # Acquiring the gene list
#  data(POA_example)
#  
#  Signature <- POA_example$POA_Rank_signature
#  
#  rowname <- get_gene_symbol(Signature)
#  
#  rownames(Signature) <- rowname$rowname
#  
#  genes <- rownames(Signature)[1:200]
#  
#  #running tisue_scMappR_custom
#  internal <- tissue_scMappR_custom(genes,Signature,output_directory = "scMappR_Test_custom", toSave = F)
#  
#  

## ----scMappR_and_pathway_analysis, eval=FALSE---------------------------------
#  
#  data(PBMC_scMappR) # load data example of PBMC bulk- and cell-sorted RNA-seq data
#  
#  bulk_DE_cors <- PBMC_example$bulk_DE_cors # 59 sex-specific DEGs in bulk PBMC (up-regulated = female-biased)
#  
#  bulk_normalized <- PBMC_example$bulk_normalized # log CPM normalized bulk RNA-seq data
#  
#  odds_ratio_in <- PBMC_example$odds_ratio_in # signature matrix developed from cell-sorted RNA-seq
#  
#  case_grep <- "_female" # flag for 'cases' (up-regulated), index is also acceptable
#  
#  control_grep <- "_male" # flag for 'control' (down-regulated), index is also acceptable
#  
#  max_proportion_change <- 10 # maximum cell-type proportion change -- this is good for cell-types that are uncomon in population and small absolute changes may yield large relative changes
#  
#  theSpecies <- "human" # these RNA-seq data have human gene symbols (and are also from human)
#  
#  # When running scMappR, it is strongly recommended to use scMappR_and_pathway analysis with the parameters below.
#  toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in,
#                                        bulk_DE_cors, case_grep = case_grep,
#                                        control_grep = control_grep, rda_path = "",
#                                        max_proportion_change = 10, print_plots = TRUE,
#                                        plot_names = "scMappR_vignette_", theSpecies = "human",
#                                        output_directory = "scMappR_vignette_",
#                                        sig_matrix_size = 3000, up_and_downregulated = TRUE,
#                                        internet = TRUE, toSave = TRUE, path = tempdir())
#  
#  

## ----tissue_ct_enrichment_example, fig.show='hide', eval=FALSE----------------
#  
#  data(POA_example)
#  POA_generes <- POA_example$POA_generes
#  POA_OR_signature <- POA_example$POA_OR_signature
#  POA_Rank_signature <- POA_example$POA_Rank_signature
#  Signature <- POA_Rank_signature
#  rowname <- get_gene_symbol(Signature)
#  rownames(Signature) <- rowname$rowname
#  genes <- rownames(Signature)[1:100]
#  
#  enriched <- tissue_by_celltype_enrichment(gene_list = genes,
#  species = "mouse",p_thresh = 0.05, isect_size = 3)
#  
#  
#  
#  

## ----process_scRNAseq_count, eval = FALSE-------------------------------------
#  
#  data(sm)
#  
#  toProcess <- list(example = sm)
#  
#  tst1 <- process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9,
#  naming_preference = "eye", rda_path = "",
#  toSave = TRUE, saveSCObject = TRUE, path = tempdir())
#  
#  
#  

## ----Seurat_Object_Generation, eval = FALSE-----------------------------------
#  
#  
#  data(sm)
#  
#  toProcess <- list(sm = sm)
#  
#  seurat_example <- process_from_count(toProcess, "test_vignette",theSpecies  = -9)
#  
#  levels(seurat_example@active.ident) <- c("Myoblast", "Neutrophil", "cardiomyoblast", "Mesothelial")
#  

## ----from_seurat_object, eval = FALSE-----------------------------------------
#  
#      generes <- seurat_to_generes(pbmc = seurat_example, test = "wilcox")
#  
#      gene_out <- generes_to_heatmap(generes, make_names = FALSE)
#  

## ----from_count_and_genes, eval = FALSE---------------------------------------
#  
#  #Create the cell-type ids and matrix
#  Cell_type_id <- seurat_example@active.ident
#  
#  count_file <- sm
#  
#  rownames_example <- get_gene_symbol(count_file)
#  
#  rownames(count_file) <- rownames_example$rowname
#  
#  # make seurat object
#  seurat_example <- process_from_count(count_file, "test_vignette",theSpecies  = "mouse")
#  
#  # Intersect column names (cell-types) with labelled CTs
#  
#  inters <- intersect(colnames(seurat_example), names(Cell_type_id))
#  
#  seurat_example_inter <- seurat_example[,inters]
#  
#  Cell_type_id_inter <- Cell_type_id[inters]
#  
#  seurat_example_inter@active.ident <- Cell_type_id_inter
#  
#  # Making signature matrices
#  
#      generes <- seurat_to_generes(pbmc = seurat_example_inter, test = "wilcox")
#  
#      gene_out <- generes_to_heatmap(generes, make_names = FALSE)
#  

## ----pathway_analysis, eval = FALSE-------------------------------------------
#  
#  data(PBMC_example)
#  bulk_DE_cors <- PBMC_example$bulk_DE_cors
#  bulk_normalized <- PBMC_example$bulk_normalized
#  odds_ratio_in <- PBMC_example$odds_ratio_in
#  case_grep <- "_female"
#  control_grep <- "_male"
#  max_proportion_change <- 10
#  print_plots <- FALSE
#  theSpecies <- "human"
#  toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in,
#                                        bulk_DE_cors, case_grep = case_grep,
#                                        control_grep = control_grep, rda_path = "",
#                                        max_proportion_change = 10, print_plots = TRUE,
#                                        plot_names = "tst1", theSpecies = "human",
#                                        output_directory = "tester",
#                                        sig_matrix_size = 3000, up_and_downregulated = FALSE,
#                                        internet = FALSE)
#  
#  twoOutFiles <- two_method_pathway_enrichment(bulk_DE_cors, "human", scMappR_vals = toOut$cellWeighted_Foldchange, background_genes = rownames(bulk_normalized), output_directory = "newfun_test",plot_names = "nonreranked_", toSave = FALSE, path=NULL)
#  
#  # The code below would save graphs and paths into the working directory.  It is commented to not make code in your working directory
#  #twoOutFiles <- two_method_pathway_enrichment(bulk_DE_cors, "human", scMappR_vals = toOut$cellWeighted_Foldchange, background_genes = rownames(bulk_normalized), output_directory = "newfun_test",plot_names = "nonreranked_", toSave = TRUE, path="./")
#  

## ----plot_barplot, eval=FALSE-------------------------------------------------
#  
#  # making an example matrix
#  term_name <- c("one", "two", "three")
#  log10 <- c(1.5, 4, 2.1)
#  
#  ordered_back_all <- as.data.frame(cbind(term_name,log10))
#  
#  #plotting
#   g <- ggplot2::ggplot(ordered_back_all, ggplot2::aes(x = stats::reorder(term_name,
#          log10), y = log10)) + ggplot2::geom_bar(stat = "identity",
#          fill = "turquoise") + ggplot2::coord_flip() + ggplot2::labs(y = "-log10(Padj)",
#          x = "Gene Ontology")
#      y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face = NULL,
#          color = "black", size = 12, angle = 35), axis.text.y = ggplot2::element_text(face = NULL,
#          color = "black", size = 12, angle = 35), axis.title = ggplot2::element_text(size = 16,
#          color = "black"))
#  
#  print(y)
#  

## ----heatmap_identification, eval=FALSE---------------------------------------
#  
#  # Generating a heatmap
#  
#  # Acquiring the gene list
#  data(POA_example)
#  
#  Signature <- POA_example$POA_Rank_signature
#  
#  rowname <- get_gene_symbol(Signature)
#  
#  rownames(Signature) <- rowname$rowname
#  
#  genes <- rownames(Signature)[1:200]
#  
#  #running tisue_scMappR_custom
#  internal <- tissue_scMappR_custom(genes,Signature,output_directory = "scMappR_Test_custom", toSave = F)
#  
#  toPlot <- internal$gene_list_heatmap$geneHeat
#  
#  
#  #Plotting the heatmap
#  
#  cex = 0.2 # size of genes
#  
#  myheatcol <- grDevices::colorRampPalette(c("lightblue", "white", "orange"))(256)
#      pheatmap::pheatmap(as.matrix(toPlot), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
#  
#  
#  

