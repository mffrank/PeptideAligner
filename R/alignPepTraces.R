
#' Aligns peptides in a pepTraces object to a reference genome
#' @importFrom pbapply pbapply
#' @return A GRangesList object with the genomic position of the found peptides

alignPepTraces <- function(pepTraces, exon, proteinseq,
                          mapIDs = FALSE, mappingTable = NULL)
{
  # Remove Decoys
  traces <- pepTraces$trace_annotation
  traces <- traces[!grepl("DECOY_", protein_id)]
  setkey(exon, tx_name)
  idmap <- exon[,.(gene_name, tx_name)]
  setkey(idmap, gene_name)
  # alignment <- lapply(traces$id, alignPeptideToGene, exon, proteinseq, pepTraces)
  # cl <- makeCluster(4)
  # alignment <- parLapply(cl,traces$id, alignPeptideToGene, exon, proteinseq, traces, idmap)
  # stopCluster(cl)
  alignment <- pbapply(traces, 1, function(trace){
    alignPeptideToGene(trace["id"], trace["protein_id"], exon, proteinseq, idmap)
  })
  names(alignment) <- traces$id
  alignment[sapply(alignment, is.null)] <- NULL
  alignment <- GRangesList(alignment)
  return(alignment)
}
