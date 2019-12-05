#' Okaty et al. (2015) Single-cell RNA-seq Data.
#'
#' We obtain the Okaty dataset [Okaty et al. (2015)](), includ- ing the raw scRNA-seq read
#' counts, the normalized gene expression values and the cell subtype labels from
#' Supplementary Material accompanying the original publication.
#'
#' @docType data
#'
#' @usage data(okaty)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @keywords datasets
#'
#' @references Okaty et al. (2015) Neuron 88 774-791
#' (\href{http://dx.doi.org/10.1016/j.neuron.2015.10.007}{DOI})
#'
#' @examples
#' data(okaty)
#' genes <- okaty$Gene.Symbol
"okaty"
