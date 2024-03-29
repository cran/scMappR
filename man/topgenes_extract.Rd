% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topgenes_extract.R
\name{topgenes_extract}
\alias{topgenes_extract}
\title{Extract Top Markers}
\usage{
topgenes_extract(generes, padj = 0.05, FC = 1.5, topNum = 30)
}
\arguments{
\item{generes}{A list of cell-type markers with fold-changes and p-values (FindMarkers output in Seurat).}

\item{padj}{The p-value (FDR) cutoff.}

\item{FC}{The fold-change cutoff.}

\item{topNum}{The number of genes to extract.}
}
\value{
\code{topgenes_extract} Returns a list of character vectors with the top (topNum) of gene markers for each cell-type. \cr
}
\description{
Internal -- Extracts strongest cell-type markers from a Seurat object.
}
\details{
Internal, this function runs through a list of outputs from FindMarkers objects in Seurat
and will extract genes past a padj and fold-change threshold. Then it extracts the topNum number of genes.
if you have not used the FindMarkers function, then a list of summary statistics with 
fold-change designated by avg_logFC and p-val by p_val_adj.
}
\examples{
 
 
 # load generes object
 data(POA_example)
 topGenes <- topgenes_extract(generes = POA_example$POA_generes)


}
