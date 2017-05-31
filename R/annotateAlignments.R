
#' Annotates protein names as the top GRangesList structure
#'
#' @param alignment Alignment output, A \code{GRangesList} with peptide sequences as names
#' @param idmap a map between peptide ID and Protein/Gene ID.
#'  Must be a \code{data.table} object with columns \code{"id", "protein_id"}.
#' @return A GRanges object with the genomic coordinates of the queried peptide
#' @import GenomicRanges

addProteinNames <- function(alignment, pep_idmap){
  stopifnot(class(alignment) == "GRangesList")
  stopifnot("data.table" %in% class(pep_idmap))
  stopifnot(c("id", "protein_id") %in% names(pep_idmap))

  setkey(pep_idmap, "protein_id")
  proteinIDs <- unique(pep_idmap$protein_id)
  res <- lapply(proteinIDs, function(x){
    pep_ids <- pep_idmap[.(x), id]
    unlist(alignment[pep_ids])
  })
  names(res) <- proteinIDs
  res <- GRangesList(res)
  return(res)
}
