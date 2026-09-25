#' Wrapper for creation, annotation and analyze of PAC object
#'
#' \code{Seqpac} Merges lanes, performs adaptor trimming and makes a first PAC object
#' with count and basic phenotypic information
#'
#' Given (at minimum) an path to input files to be analyzed, PAC_create will wrap around
#' all necessary preparatory steps (\code{\link{merge_lanes}}, \code{\link{make_counts}}, and
#' \code{\link{make_PAC}}) to produce a PAC object. This wrapper comes with the caveat that
#' it uses default values in most functions. If that does not suit your analysis, we recommend
#' you to run each function on its own!
#'
#' @family PAC generation
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#' @param lanes allows to define whether lanes should be merged with \code{\link{merge_lanes}}
#'   (may be necessary after using e.g BaseSpace to download fasta files from Illumina cloud).
#'   Default=NULL, where no lane merging will be performed.
#'
#' @param trim Character vector to define whether trimming of adaptors should be performed.
#'   Current accepted values are "default_neb" and "default_illumina", trimming
#'   adaptors from library preparation with NebNext or Illumina technologies
#'   respectively. More advanced and manual adaptor trimming is possible with
#'   \code{\link{make_trim}}. Default=NULL, where no adaptor trimming will be performed.
#'
#' @param input Character indicating path to where input files can be found in fasta 
#'   or fastq format. Depending on value of lanes and trim, the function assumes whether
#'   to expect files in path to be trimmed or not and whether to merge lanes.
#'   
#' @param output Character indicating path to where output files will be stored, 
#'   if applicable (such as in merge_lanes). Default=NULL. 
#'   
#' @param pheno Can be either a data frame, or a string to directory for a comma
#'  separated .csv file that will be used to produce the PAC object. Default=NULL,
#'  where progress report from counts will be added to pheno. 
#'   
#' @param input_genome Character indicating path to reference genome in fasta (.fa) 
#'  format to use for annotation by \code{\link{map_reanno}} with import="genome". 
#'  Reference genome should have a bowtie index, please see 
#'  ??Rbowtie::bowtie_build for instructions. Default=NULL, where no genome-based
#'  mapping will be performed.
#'
#' @param input_biotype Character indicating path to reference fasta (.fa) file for 
#'  bioinformatic annotation by \code{\link{map_reanno}} with import="biotype". 
#'  This wrapper is justerad to the ncRNA fasta from Ensembl, with biotypes such 
#'  as rRNA, tRNA and miRNA available. Reference fasta should have a bowtie index, 
#'  please see ??Rbowtie::bowtie_build for instructions. Within this wrapper, a 
#'  hierarcy is made with \code{\link{simplify_reanno}}, with the following 
#'  hiearchy: rRNA, miRNA, tRNA, snoRNA, snRNA. If other references wish to be 
#'  used, please run each function separately. Default=NULL, where no biotype-
#'  based mapping will be performed.
#' 
#' @param override Logical whether or not the function should prompt you for a
#'  question if there are files in output. As default, override=FALSE will
#'  prevent deleting large files by accident, but requires an interactive R
#'  session. Setting override=TRUE may solve non-interactive problems. 
#'   
#' @param pheno_target List with: 1st object being a character vector
#'  of target column in Pheno and 2nd object being a character vector of the target
#'  group(s) in the target Pheno column (1st object).
#' 
#' @param norm Character indicating what type of data to be used. If 'counts',
#'  the raw counts in Counts will be used (default). If "cpm", cpm values will 
#'  be used. As of now, only "counts" and "cpm" may be used for this wrapper.
#' 
#' @param anno_target Character vector with the name of the target column in
#'  Anno or the name of the annotation column in case of input being a
#'  dataframe.
#'  
#' @param filter_hit If filter_hit is TRUE, sequences without a hit to the reference
#'  genome(s) defined in input_genome up until the mismatches defined will be removed.
#'  Default is TRUE.
#' 
#' @param model Character of model used to run \code{PAC_deseq}.
#' 
#' @param ... Arguments to be passed on to \code{\link{PAC_create}}, \code{\link{PAC_map}},
#'  and \code{\link{PAC_annotate}}.
#'   
#' @return PAC object
#'   
#' @examples
#' 
#' ###########################################################
#' ##----------------------------------------
#' 
#' #load in the test reference files and ensure correct Bowtie files are available
#' 
#' input = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' 
#' ## tRNA:
#' trna_file <- system.file("extdata/trna", "tRNA.fa",
#'                          package = "seqpac", mustWork = TRUE)
#' trna_dir <- dirname(trna_file)
#' 
#' if(!sum(stringr::str_count(list.files(trna_dir), ".ebwt")) ==6){
#'   Rbowtie::bowtie_build(trna_file,
#'                         outdir=trna_dir,
#'                         prefix="tRNA", force=TRUE)
#' }
#' ## rRNA:
#' rrna_file <- system.file("extdata/rrna", "rRNA.fa",
#'                          package = "seqpac", mustWork = TRUE)
#' rrna_dir <- dirname(rrna_file)
#' 
#' if(!sum(stringr::str_count(list.files(rrna_dir), ".ebwt")) ==6){
#'   Rbowtie::bowtie_build(rrna_file,
#'                         outdir=rrna_dir,
#'                         prefix="rRNA", force=TRUE)
#' }
#' 
#' input_biotype <- list(trna= trna_file, rrna= rrna_file)
#' 
#' ## "Genome" (here we reuse the rRNA reference just as an example):
#' input_genome <- list(genome=rrna_file)
#' 
#' 
#' #Manually making a pheno for this example
#' pheno <- data.frame(row.names=(c("seqpac_fq1", "seqpac_fq2",
#'                                  "seqpac_fq3", "seqpac_fq4", "seqpac_fq5", "seqpac_fq6")),
#'                     stage=rep(c("Stage1", "Stage3"), each=3))
#' 
#' Seqpac(input=input, input_genome=input_genome, input_biotype=input_biotype, 
#'        pheno_target = list("stage"), anno_target=list("Biotypes"), 
#'        pheno=pheno, norm="cpm", filter_hit=FALSE, override=TRUE)
#'
#' @export


Seqpac <- function(lanes=NULL, trim=NULL, input, output=NULL, input_genome=NULL, 
                   input_biotype=NULL, pheno_target=NULL, norm=NULL, 
                   anno_target=NULL, model=NULL, override=TRUE,
                   filter_hit=FALSE, 
                   pheno=NULL, ...)
  {
  cat("Creating PAC object ... \n")
  pac <- PAC_create(lanes=lanes, trim=trim, input=input, output=output, pheno=pheno, ...)
  cat("PAC created. \n")
  print(pac)
  cat("Annotating PAC object ... \n")
  pac <- PAC_map(input_genome=input_genome, input_biotype=input_biotype, output=output, 
                 PAC=pac, override=override, filter_hit=filter_hit, ...)
  cat("Analyzing PAC object ... \n")
  if(!is.null(norm)){
    cat("Normalizing PAC object with:", print(norm))
    pac <- PAC_norm(pac, norm=norm)
  }
  result_list <- PAC_analyze(PAC=pac, pheno_target=pheno_target, norm=norm,
                     anno_target=anno_target, model=model, output=output)
  closeAllConnections()
  return(pac)
}
