#' Aligns a provided peptide sequence to a provided Gene sequence.
#'
#' @param peptide Peptide sequence as character
#' @param geneID The geneID from which the peptide originated
#' @param exon A \code{data.table} with exon information, obtained by a download function(see ...)
#' @param proteinseq a table with protein sequences for every transcript in the reference genome.
#' @param idmap a map between geneID and Transcript ID. Must be a keyed \code{data.table} object.
#' @return A GRanges object with the genomic coordinates of the queried peptide
#' @examples
#' alignPeptideToGene(peptide = "LHAVVETLVNHR", geneID = "ENSG00000028528",
#'                    exon = exon_dt, proteinseq, idmap)
#' @import data.table
#' @export
alignPeptideToGene <- function(peptide,geneID, exon, proteinseq, idmap){
  # geneID <- subset(traces, id == peptide)$protein_id
  # transcriptIDs <- subset(exon, gene_name == geneID)$tx_name
  transcriptIDs <- idmap[.(geneID), tx_name]
  if(!is.na(transcriptIDs[1])){

    for(transcript in transcriptIDs){
      peptide_loc <- alignPeptideToTranscript(peptide, transcript, exon, proteinseq)
      if(!is.null(peptide_loc)){
        return(peptide_loc)
      }
    }
    return(NULL)
  }else{
    #@TODO implement genomewide alignment
    return(NULL)
  }
}




#' Aligns a provided peptide sequence to a provided Transcript.
#'
#' @param peptide Peptide sequence as character
#' @param transcriptID The transcriptID from which the peptide possibly originated
#' @param exon A \code{data.table} with exon information, obtained by a download function(see ...)
#' @param proteinseq a table with protein sequences for every transcript in the reference genome.
#' @return A GRanges object with the genomic coordinates of the queried peptide

#' @importFrom  Biostrings matchPattern
#' @import GenomicRanges

alignPeptideToTranscript <- function(peptide, transcriptID, exon, proteinseq){

  # peptide <- pepTraces$trace_annotation$id[a]
  peptide_seq <- gsub("\\(.*?\\)", "", peptide)
  protein <- subset(proteinseq, tx_name == transcriptID)$peptide
  if(!is.na(protein[1])){
    m <- matchPattern(peptide_seq, AAString(protein))
  }else{
    return(NULL)
  }

  if(!is.na(start(m)[1])){
    peptide_loc <- GRanges()

    for(i in 1:length(m)){
      pep_cds_start <- 3*start(m)[i]-2
      pep_cds_end <- 3*end(m)[i]
      #Binary search for increased performance
      tr_exons <- exon[.(transcriptID)]
      done = F
      while(done == F){
        # ex_one <- exon[which(exon$tx_name == transcriptID & exon$cds_start <= cds_start & exon$cds_end >= cds_start),]
        ex_one <- subset(tr_exons,(cds_start <= pep_cds_start) & (cds_end >= pep_cds_start))
        chr_start <- ex_one$cds_chr_start + pep_cds_start - ex_one$cds_start
        if(ex_one$cds_end >= pep_cds_end){
          chr_end <- ex_one$cds_chr_start + pep_cds_end - ex_one$cds_start
          peptide_loc <- append(peptide_loc, GRanges(ex_one$chromosome_name, strand = ex_one$strand,
                                                     IRanges(start = chr_start, end = chr_end)))
          done = T
        }else{
          chr_end <- ex_one$cds_chr_end
          peptide_loc <- append(peptide_loc, GRanges(ex_one$chromosome_name, strand = ex_one$strand,
                                                     IRanges(start = chr_start, end = chr_end)))
          pep_cds_start <- ex_one$cds_end + 1
        }
      }
    }
    return(peptide_loc)
  }else{
    return(NULL)
  }
}
