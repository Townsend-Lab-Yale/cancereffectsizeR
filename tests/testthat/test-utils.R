test_that('MNV detection and identify_maf_variants()', {
  
  # Check MNV detection
  complex_variants = data.table(Chromosome = c('9', '9', '12', '12', '14', '14', '1', '1'), 
                                Start_Position = c(3.5e4, 3.5e4, 1e6, 1e6, 5e7, 5e7 + 2, 1e6, 1e6 + 2), 
                                Reference_Allele = c('C', 'C', 'G', 'G', 'T', '-', 'T', 'G'), 
                                Tumor_Allele = c('G', 'T', 'T', 'A', 'G', 'AA', 'C', 'A'), 
                                Unique_Patient_Identifier = c('a', 'a', 'b', 'b', 'c', 'c', 'd', 'd'))
  preloaded = preload_maf(complex_variants, refset = 'ces.refset.hg19')
  expect_equal(preloaded[! is.na(problem), unique(problem)], 'merged_with_nearby_variant')
  preloaded = preloaded[is.na(problem)]
  
  expected = structure(list(Chromosome = c("9", "12", "14", "1"),
                            Start_Position = c(35000, 1e+06, 5e+07, 1e+06),
                            Reference_Allele = c("C,C", "G,G", "T,(+2)-", "T,(+2)G"),
                            Tumor_Allele = c("G,T", "T,A", "G,AA", "C,A"),
                            Unique_Patient_Identifier = c("a", "b", "c", "d"),
                            variant_type = c("other", "other", "other", "other"),
                            variant_id = c("9:35000_C>G,9:35000_C>T", "12:1000000_G>T,12:1000000_G>A", 
                                           "14:50000000_T>G,14:50000002_ins_AA", "1:1000000_T>C,1:1000002_G>A"
                            )))
  expect_equal(preloaded[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele, Unique_Patient_Identifier,
                             variant_type, variant_id)],
               expected, check.attributes = F)
  reidentified = identify_maf_variants(copy(preloaded[, -"variant_id"]))
  expect_equal(reidentified[, .(variant_id, variant_type)], preloaded[, .(variant_id, variant_type)])
  
})

test_that('identify_maf_variants marks illegal on complex variant with unbalanced ref/alt', {
          bad_complex = data.table(Chromosome = '1', Start_Position = 5, Reference_Allele = c('A,C'), 
                                   Tumor_Allele = c('G,T,C'),
                                   Unique_Patient_Identifier = 'a')
          expect_equal(identify_maf_variants(bad_complex)[, unlist(.(variant_type, variant_id))],
                    c("illegal", NA_character_))
          expect_equal(preload_maf(bad_complex, refset = 'ces.refset.hg19')[, problem], 'invalid_record')
})

test_that('identify_maf_variants does not mangle repeated complex identifications,', 
          {
            input_maf = data.table(Unique_Patient_Identifier = c("patient", "patient", "patient"), 
                                   Chromosome = c("chr8", "chr8", "chr8"), Start_Position = 96287061:96287063, 
                                   Reference_Allele = c("T", "C", "T"), Tumor_Allele = c("C", "A", "C"))
            complex_maf = preload_maf(maf = input_maf, refset = "ces.refset.hg38")
            expect_equivalent(complex_maf[is.na(problem), unlist(.(Unique_Patient_Identifier, Chromosome, Start_Position, 
                                                       Reference_Allele, Tumor_Allele, variant_id))],
                              c("patient", "8", "96287061", "T,(+1)C,(+2)T", "C,A,C", 
                                "8:96287061_T>C,8:96287062_C>A,8:96287063_T>C"))
            reidentified = identify_maf_variants(copy(complex_maf))
            expect_equivalent(reidentified, complex_maf)
          })




