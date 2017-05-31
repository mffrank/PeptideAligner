
#' Annotates protein names as the top GRangesList structure
#'
#' @param alignment Alignment output, A \code{GRangesList} with peptide sequences as names
#' @param idmap a map between peptide ID and Protein/Gene ID.
#'  Must be a \code{data.table} object with columns \code{"id", "protein_id"}.
#' @return A GRanges object with the genomic coordinates of the queried peptide
#' @import GenomicRanges

addProteinNames <- function(alignment, idmap){
  stopifnot(class(alignment) == "GRangesList")
  stopifnot(class(idmap) %in% "data.table")
  stopifnot(c("id", "protein_id") %in% names(idmap))

  res <- GRangesList()
  setkey(idmap, "protein_id")
  res <- lapply(idmap$protein_id, function(x){
    pep_ids <- idmap[.(x), id]
    alignment[pep_ids]
  })
  res <- do.call("append", res)
  return(res)
}
