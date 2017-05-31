# Whole outout data from initial analysis
data_out <- readRDS("data-raw/pepAlignment_customDB_semitryptic_rev_decoys.rda")


#Input data (too large to store here)
library(data.table)
annotation_path <- "Y:/Master_Project/data/Human_genome/GRCh38/CustomProDB_annotation"
load(paste(annotation_path, "/exon_anno.RData", sep = ""))
load(paste(annotation_path,"/proseq.RData", sep = ""))
pepTraces <- readRDS(("Y:/Master_Project/data/HEK293/Peptide_Traces/customDB_semitryptic_rev_decoys/pepTraces_Iso_collapsed.rda"))

# Arbitrary peptide for testing
pep1 <- "YSVLLPTYNER"
prot1 <- pepTraces$trace_annotation[id == pep1, protein_id]
pep2 <- "TVFQFDVQR"
prot2 <- pepTraces$trace_annotation[id == pep2, protein_id]
testpeps <- pepTraces$trace_annotation[protein_id %in% c(prot1, prot2), id]
testtranscripts <- exon[exon$gene_name %in% c(prot1, prot2), "tx_name"]

pepTraces_test <- subset(pepTraces, testpeps)

exon_test <- as.data.table(exon)
exon_test <- exon_test[gene_name %in% c(prot1, prot2)]

proteinseq_test <- as.data.table(proteinseq)
proteinseq_test <- proteinseq_test[tx_name %in% testtranscripts]

devtools::use_data(pepTraces_test, proteinseq_test, exon_test)
