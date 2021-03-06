#' PBMC_scMappR
#' 
#' Toy example of data where cell-weighted fold-changes and related downsteam analyses can be completed.
#'
#' A named list called "PBMC_example" containing the count data, signature matrix, and DEGs. The count data and signature matrix are shortened to fit the size of the package and do not reflect biologically relevant data.
#'
#' @rdname PBMC_example
#' @name PBMC_example
#'
#' @usage data(PBMC_example)
#'
#' @format A list containing three data frames, normalized count data, a signature matrix, and a list of differentially expressed genes.
#' \describe{
#'   \item{bulk_normalized}{ A 3231 x 9 matrix where rows are genes, columns are samples, and the matrix is filled with CPM normalized counts.}
#'   \item{odds_ratio_in}{ A 2336 x 7 matrix where rows are genes, columns are cell-types and matrix is filled with the odds-ratio that a gene is in each cell-type.}
#'   \item{bulk_DE_cors}{A 59 x 3 matrix of sex-specific genes found between male and female PBMC samples (female biased = upregulated). row and rownames are genes, columns are gene name, FDR adjusted p-value, and log2 fold-change. DEGs were computed with DESeq2 and genes with a log2FC > 1 were kept. }
#' }
#' @examples 
#' data(PBMC_example)
#' @keywords datasets
#' 
NULL
