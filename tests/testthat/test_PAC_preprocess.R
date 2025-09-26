context("PAC preprocessing and summary\n")
library(seqpac)

invisible(capture.output(
    test_that("Testing PAC_filter and  PAC_filtsep ...", {
      load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                       package = "seqpac", mustWork = TRUE))
      
      pac_lowfilt <- PAC_filter(pac, nucleotide_range=c(15,80), threshold=5, 
                                coverage=20, norm = "counts",  
                                pheno_target=NULL, anno_target=NULL)
      
      pac_subset <- PAC_filter(pac_lowfilt, subset_only = TRUE,
                               pheno_target=list("batch", c("Batch1", "Batch2")), 
                               anno_target=list("Size", "22"))
      
      filtsep <- PAC_filtsep(pac_lowfilt, norm="counts", threshold=, 
                             coverage=100, pheno_target= list("stage"),
                             output="binary")
      
      pac_filt <- PAC_filter(pac_lowfilt, subset_only = TRUE,
                             anno_target= row.names(filtsep))
      
      expect_true(PAC_check(pac_lowfilt))
      expect_true(PAC_check(pac_subset))
      expect_true(PAC_check(pac_filt))
      
    })))

invisible(capture.output(    
    test_that("Testing PAC_norm and PAC_summary ...", {
      load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                       package = "seqpac", mustWork = TRUE))
      pac  <- PAC_norm(pac, norm="cpm")
      suppressMessages(   
          pac  <- PAC_norm(pac, norm="vst", merge_pac = TRUE)
      )
      suppressMessages( 
          pac  <- PAC_norm(pac, norm="rlog")
      )
      
      df_norm  <- PAC_norm(pac, norm="cpm", merge_pac = FALSE)
      
      pac  <- PAC_summary(pac, norm="cpm", pheno_target=list("stage"))
      pac  <- PAC_summary(pac, norm="cpm", pheno_target=list("batch", c("Batch1", "Batch2")))
      df_sum  <- PAC_summary(pac, norm="cpm", pheno_target=list("batch", c("Batch1", "Batch2")), merge_pac = FALSE)
     
      expect_true(PAC_check(pac))
      expect_equal(norm(pac)$cpm, df_norm)
      expect_equal(summary(pac)$cpmMeans_batch, df_sum[[1]])
      
    })))
  



