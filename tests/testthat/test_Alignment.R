context("Alignment")

test_that("Alignment Output of transcript",{
  al <- alignPeptideToGene(peptide = pepTraces_test$trace_annotation$id[1],
                     geneID = pepTraces_test$trace_annotation$protein_id[1],
                     exon = exon_test, proteinseq = proteinseq_test,
                     idmap = idmap_test)

  expect_equal(start(al), c(50940869, 50936148))
  expect_equal(end(al), c(50940933, 50936148))
})
