################################################################################
###-----------------------------------------------------------------------------
### This file contains all new S4 classes, generics and methods for seqpac. This
### involves:
### 1. The PAC object 
### 2. The reanno object 
### Generating these objects are controlled by (1) make_PAC and (2) make_reanno
### functions. From these functions, S3 versions of the objects can be
### generated using the 'output=' option. In addition S3 versions are available
### through coercion methods, eg. as(PAC,"S3"). The map object generated from
### the PAC_mapper function is a regular list.
################################################################################
###-----------------------------------------------------------------------------
### the S4 PAC object (make_PAC):
##
#' S4 object for PAC 
#' 
#' S4 version of the PAC object. 
#' 
#' Just like the S3 version, all sub-tables must have identical sequence (row)
#' names. To extract tables from S4 objects, use the getter commands or use @
#' (e.g. PAC@Pheno) just like $ for the S3 version (e.g. pheno(PAC)). You may
#' also use a coercion method, e.g. \code{PAC_S3 <- as(PAC_S4, "list")}
#' 
#' Holds up to 5 slots: \describe{
#'    \item{\strong{Pheno}}{data.frame > information about the samples}
#'    \item{\strong{Anno}}{data.frame > sequences annotations}
#'    \item{\strong{Counts}}{data.frame > sequence counts over samples}
#'    \item{\emph{norm (optional)}}{list of data.frames > normalized counts}
#'    \item{\emph{summary (optional)}}{
#'    list of data.frames > summarized counts or normalized counts}
#'    }
#'    
#' @rdname PAC
#' @return Generates a S4 PAC-object.
#' @examples
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'                   
#' # check type
#' class(pac)
#' isS4(pac)
#'  
#'  
#' # Turns S4 PAC object into a S3
#' pac_s3 <- as(pac, "list")
#' 
#' # Turns S3 PAC object into a S4
#' pac_s4 <- as.PAC(pac_s3)
#' 
#' 
#' @export
#' 
setClass(Class = "PAC",
          slots=c(Pheno = "data.frame", 
                  Anno= "data.frame",
                  Counts= "data.frame",
                  norm= "list",
                  summary="list"
                  ))

#-------------------------------
# Constructor PAC
#' @rdname PAC
#' @param Pheno Phenotype data.frame table with sample ID as row names and where
#'   columns are variables associated with the samples. Can be generated using
#'   the make_pheno function.
#' @param Anno Annotation data.frame table with unique sequences as row names
#'   and where columns are variables associated with the sequences.
#' @param Counts Counts table (data.frame) with unique sequences as row names
#'   (identical to Anno) and where columns are sample ID (identical to row names
#'   for Pheno. Should contain raw counts and can be generated using the
#'   make_counts function.
#' @param norm List of data.frames. May be regarded as a "folder" with
#'   normalized counts tables. The listed data.frames must have identical
#'   sequence names as Counts/Anno (rows) and sample IDs (columns) as
#'   Counts/Pheno. Can be generated using the PAC_norm function, but seqpac
#'   functions will attempt use any normalized table stored in norm that are in
#'   agreement with the above cafeterias.
#' @param summary List of data.frames, just like norm, but contains summarized
#'   raw or normalized counts (e.g. means, standard errors, fold changes).
#'   Important, the listed data.frames must have identical sequence names as
#'   Counts/Anno (rows), but columns don't need sample IDs. Can be generated
#'   using the PAC_summary function, but seqpac functions will attempt use any
#'   summarized table stored in the summary "folder".
#' @export
PAC <- function(Pheno, Anno, Counts, norm, summary){
  new("PAC", Pheno=Pheno, Anno=Anno, Counts=Counts,
      norm=norm, summary=summary)}

#-------------------------------
# Test validity
setValidity("PAC", function(object){
   # Is it empty
   if(is.null(length(object))){
     return("Empty object")
   }
   # Test 1
   pac <- as(object, "list")
   return(PAC_check(pac))
 })


#-------------------------------
# Coercion method PAC
setAs("PAC", "list",
      function(from){
        pac <- list(Pheno=seqpac::pheno(from),
                    Anno=seqpac::anno(from),
                    Counts=seqpac::counts(from))
        
        if(!is.null(seqpac::norm(from)[[1]])){
         pac <- c(pac, list(norm= seqpac::norm(from)))
        }
        if(!is.null(seqpac::summary(from)[[1]])){
         pac <- c(pac, list(summary= seqpac::summary(from)))
        }
        class(pac) <- c("PAC_S3","list")
        return(pac)
        })

###############################################
#' Converts an S3 PAC into a S4 PAC
#'
#' \code{as.PAC} Converts an S3 PAC object (list) to an S4 PAC object
#' 
#' Seqpac comes with two versions of the PAC object, either S4 or S3. The S3 PAC
#' is simply a list where each object can be received using the $-sign. The S4
#' version is a newer type of R object with predefined slots, that is received
#' using the @-sign. This function converts an S3 PAC object (list) to an S4 PAC
#' object. You can also use the S4 coercion method "as" to turn an S4 PAC into
#' an S3. See examples below.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/Danis102} for updates on the current
#'   package.
#'
#' @param from S3 PAC-object containing at least a Pheno table with samples as
#'   row names, Anno table with sequences as row names and a Count table with
#'   raw counts (columns=samples, rows=sequences).
#' 
#' @return An S4 PAC object.
#'   
#' @examples
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#'                   
#' # check type
#' class(pac)
#' isS4(pac)
#'  
#'  
#' # Turns S4 PAC object into a S3
#' pac_s3 <- as(pac, "list")
#' 
#' # Turns S3 PAC object into a S4
#' pac_s4 <- as.PAC(pac_s3)
#' 
#' 
#' @export
#' 
as.PAC <- function(from){
       pac <- seqpac::PAC(Pheno=from$Pheno,
                   Anno=from$Anno,
                   Counts=from$Counts,
                   norm=list(NULL),
                   summary=list(NULL)
                   )
       if("norm" %in% names(from)){
              pac@norm <- from$norm
        }
       if("summary" %in% names(from)){
              pac@summary <- from$summary
        }
       return(pac)
        }


###############################################################################
###############################################################################
###############################################################################
###----------------------------------------------------------------------------
### the S4 reanno object (make_reanno):
#'
#' S4 object for reanno output 
#' 
#' Holds the imported information from mapping using the map_reanno function.
#' All information are held in tibble class (tbl/tbl_df) tables from the tibble
#' package. 
#' 
#' @return Contains two slots:  1. Overview: A table holding a summary of the
#' mapping. 2. Full_anno: Lists of tables holding the full imports per mismatch
#' cycle (mis0, mis1 etc) and reference (mi0-ref1, mis0-ref2, mis1-ref1,
#' mis1-ref2 etc).
#' 
#' @rdname reanno
#' 
#' @importClassesFrom tibble tbl_df
#' 
#' @examples
#' 
#' ##### Create an reanno object
#' 
#' ##  First load a PAC- object
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#'  
#' ##  Then specify paths to fasta references
#' # If you are having problem see the vignette small RNA guide for more info.
#'  
#'  trna_file <- system.file("extdata/trna", "tRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)
#'  trna_dir <- dirname(trna_file)
#'  
#'  if(!sum(stringr::str_count(list.files(trna_dir), ".ebwt")) ==6){     
#'      Rbowtie::bowtie_build(trna_file, 
#'                            outdir=trna_dir, 
#'                            prefix="tRNA", force=TRUE)
#'                            }
#'  
#'  ref_paths <- list(trna= trna_file)
#'                                     
#' ##  Add output path of your choice.
#' # Here we use the R temporary folder depending on platform
#'  output <- paste0(tempdir(),"/seqpac/test")
#' 
#' ## Make sure it is empty (otherwise you will be prompted for a question)
#' out_fls  <- list.files(output, recursive=TRUE)
#' closeAllConnections()
#' suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
#'
#' ##  Then map your PAC-object against the fasta references
#'  map_reanno(pac, input=ref_paths, output=output,
#'                type="internal", mismatches=0,  import="biotype", 
#'                threads=2, keep_temp=FALSE, override=TRUE)
#' 
#' ##  Then import and generate a reanno-object of the temporary bowtie-files
#' reanno <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
#'                                     
#' ## Now make some search terms against reference names to create shorter names
#' # Theses can be used to create factors in downstream analysis
#' # One search hit (regular expressions) gives one new short name 
#' bio_search <- list(trna =c("_tRNA", "mt:tRNA"))
#'  
#'  
#' ## You can merge directly with your PAC-object by adding your original 
#' # PAC-object, that you used with map_reanno, to merge_pac option.
#'  pac <- add_reanno(reanno, bio_search=bio_search, 
#'                        type="biotype", bio_perfect=FALSE, 
#'                        mismatches = 0, merge_pac=pac)
#'                        
#'                        
#' ## Turn your S3 list to an S4 reanno-object
#' isS4(reanno)
#' names(reanno)
#' table(overview(reanno)$trna)
#'  
#' reanno_s3 <- as(reanno, "list")
#' isS4(reanno_s3)   
#'
#' # Turns S3 reanno object into a S4 
#' reanno_s4 <- as.reanno(reanno_s3)
#' isS4(reanno_s4) 
#'  
#' @export

setClass(Class = "reanno",
          slots=c(Overview = "tbl_df", 
                  Full_anno= "list"))

#-------------------------------
# Constructor reanno
#' @rdname reanno
#' @param Overview A tibble data.frame with summarized results from mapping
#'   using the \code{\link{map_reanno}} function that has been imported into R
#'   using the \code{\link{make_reanno}} function. Rows represents sequences
#'   from the original PAC-object.
#' @param Full_anno A multi-level list with tibble data.frames that contains all
#'   that was imported by \code{\link{make_reanno}}
#' @export
reanno <- function(Overview, Full_anno){
  new("reanno", Overview=Overview, Full_anno=Full_anno)}

#-------------------------------
# Test validity reanno

setValidity("reanno", function(object){
  # Is it empty
  if(is.null(length(object))){
    return("Empty object")
  }
  # Test 1
  seq_nam <- dplyr::pull(overview(object)[1])
  test1 <- lapply(full(object), function(y){
    lapply(y, function(z){identical(seq_nam, dplyr::pull(z[1]))})
  })
  logi_test1   <- unlist(test1)
  # Test 2
  test2 <- lapply(full(object), function(y){
    lapply(y, function(z){methods::is(z, c("tbl_df"))})
    })
  logi_test2   <- unlist(test2)
  # Return
  if(any(!logi_test1)){ 
    return(cat("\nSequence names are not identical.",
               "\nPlease check @Full_anno tables:\n",
               names(logi_test1)[which(!logi_test1)]))
  }
  if(any(!logi_test2)){ 
    return(cat("\n@Full_anno is not a tibble.",
               "\nPlease check @Full_anno tables:\n",
               names(logi_test2)[which(!logi_test2)]))
  }
  TRUE
})

#-------------------------------
# Coercion method reanno
#
setAs("reanno", "list",
      function(from){
        rn <- list(Overview=overview(from), 
                      Full_anno=full(from))
        class(rn) <- c("reanno_S3", "list")
        return(rn)
        })


###############################################
#' Converts an S3 reanno into a S4 reanno
#'
#' \code{as.reanno} Converts an S3 reanno object (list) to an S4 reanno object
#' 
#' Seqpac comes with two versions of the reanno object, either S4 or S3. The S3 
#' is simply a list where each object can be received using the $-sign. The S4
#' version is a newer type of R object with predefined slots, that is received
#' using the @-sign. This function converts an S3 reanno object (list) to an S4 
#' object. You can also use the S4 coercion method "as" to turn an S4 into
#' an S3. See examples below.
#' 
#' @family PAC reannotation
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param from S3 reanno object.
#' 
#' @return An S4 reanno object.
#'   
#' @examples
#' 
#' ##### Create an reanno object
#' 
#' ##  First load a PAC- object
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#'  
#' ##  Then specify paths to fasta references
#' # If you are having problem see the vignette small RNA guide for more info.
#'  
#'  trna_file <- system.file("extdata/trna", "tRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE)
#'  trna_dir <- dirname(trna_file)
#'  
#'  if(!sum(stringr::str_count(list.files(trna_dir), ".ebwt")) ==6){     
#'      Rbowtie::bowtie_build(trna_file, 
#'                            outdir=trna_dir, 
#'                            prefix="tRNA", force=TRUE)
#'                            }
#'  
#'  ref_paths <- list(trna= trna_file)
#'                                     
#' ##  Add output path of your choice.
#' # Here we use the R temporary folder depending on platform
#'  output <- paste0(tempdir(),"/seqpac/test")
#' 
#' ## Make sure it is empty (otherwise you will be prompted for a question)
#' out_fls  <- list.files(output, recursive=TRUE)
#' closeAllConnections()
#' suppressWarnings(file.remove(paste(output, out_fls, sep="/")))
#'
#' ##  Then map your PAC-object against the fasta references
#'  map_reanno(pac, input=ref_paths, output=output,
#'                type="internal", mismatches=0,  import="biotype", 
#'                threads=2, keep_temp=FALSE, override=TRUE)
#' 
#' ##  Then import and generate a reanno-object of the temporary bowtie-files
#' reanno <- make_reanno(output, PAC=pac, mis_fasta_check = TRUE)
#'                                     
#' ## Now make some search terms against reference names to create shorter names
#' # Theses can be used to create factors in downstream analysis
#' # One search hit (regular expressions) gives one new short name 
#' bio_search <- list(trna =c("_tRNA", "mt:tRNA"))
#'  
#'  
#' ## You can merge directly with your PAC-object by adding your original 
#' # PAC-object, that you used with map_reanno, to merge_pac option.
#'  pac <- add_reanno(reanno, bio_search=bio_search, 
#'                        type="biotype", bio_perfect=FALSE, 
#'                        mismatches = 0, merge_pac=pac)
#'                        
#'                        
#' ## Turn your S3 list to an S4 reanno-object
#' isS4(reanno)
#' names(reanno)
#' table(overview(reanno)$trna)
#'  
#' reanno_s3 <- as(reanno, "list")
#' isS4(reanno_s3)   
#'
#' # Turns S3 reanno object into a S4 
#' reanno_s4 <- as.reanno(reanno_s3)
#' isS4(reanno_s4) 
#'  
#' @export
#' 
as.reanno <- function(from){
       reanno <- seqpac::reanno(Overview=from$Overview,
                   Full_anno=from$Full_anno
                   )
       return(reanno)
        }



################################################################################
###-----------------------------------------------------------------------------




