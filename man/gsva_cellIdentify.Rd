% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsva_cellIdentify.R
\name{gsva_cellIdentify}
\alias{gsva_cellIdentify}
\title{Cell-type naming with GSVA}
\usage{
gsva_cellIdentify(
  pbmc,
  theSpecies,
  naming_preference = -9,
  rda_path = "",
  toSave = FALSE
)
}
\arguments{
\item{pbmc}{Processed Seurat object without named cells.}

\item{theSpecies}{"human" or "mouse" -- it will determine which species cell-type markers will originate from.}

\item{naming_preference}{Once top cell-type markers are identified, naming_preferences will then extract CT markers within a more appropriate tissue type.}

\item{rda_path}{Path to pre-computed cell-type .gmt files (rda objects).}

\item{toSave}{If scMappR is allowed to write files and directories.}
}
\value{
List with the following elements:
\item{cellMarker}{Most likely cell-types predicted from CellMarker database.}
\item{panglao}{Most likely cell-types predicted from Panglao database.}
\item{avg_expression}{Average expression of each gene in each cell-type.}
}
\description{
This function computes the mean expression of every cell-type before predicting the most likely cell-type using the GSVA method.
}
\details{
This function inputs a Seurat object and uses the average normalized expression of each gene in each cluster to identify cell-types using the GSVA method.
}
\examples{
\donttest{

data(sm)
toProcess <- list(example = sm)
tst1 <- process_from_count(countmat_list = toProcess,name = "testProcess",
                            theSpecies = "mouse")
cellnames <- gsva_cellIdentify(pbmc = tst1, theSpecies = "mouse",
 naming_preference = "brain", rda_path = "")

}
 
}
