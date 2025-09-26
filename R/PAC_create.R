#' Wrapper for creation of PAC object
#'
#' \code{PAC_create} Merges lanes, performs adaptor trimming and makes a first PAC object
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
#'
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
#' @return PAC object with values in pheno and count
#'   
#' @examples
#' 
#' ###########################################################
#' ##----------------------------------------
#' 
#' # First create an annotation blank PAC with group means
#' sys_path = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' fq <- list.files(path = sys_path, pattern = "fastq", all.files = FALSE,
#'                 full.names = TRUE)
#'  
#' # Create an output folder
#' input <- paste0(tempdir(), "/lanes/")
#' output <- paste0(tempdir(), "/merged/")
#' dir.create(input, showWarnings=FALSE)
#' dir.create(output, showWarnings=FALSE)
#'
#' # Fix compatible file names
#' file.copy(from = fq, to = input)
#' old_fls <- list.files(input, full.names=TRUE)
#' new_sample <- c(rep("sample1_", times=3), rep("sample2_", times=3))
#' new_lane <- rep(c("lane1","lane2","lane3"), times=2)
#' new_fls <- paste0(input,new_sample, new_lane, ".fastq.gz")
#' file.rename(from = old_fls, to = new_fls)
#' 
#' # Then merge the fastq files
#' pac <- PAC_create(lanes=TRUE, trim="default_neb", input, output)
#'
#'
#' @export



PAC_create <- function(lanes=NULL, trim=NULL, input, output=NULL){
  
  inpath<-input
  outpath<-output
  
  if(!is.null(lanes)){
    cat("Will now merge lanes...")
    merge_lanes(input=inpath, output = outpath)
  }
  if(!is.null(lanes)){
    
    if(!is.null(trim)){
      cat("Will now merge lanes...")
      if(trim=="default_neb"){
        counts<-make_counts(input=outpath,
                            trimming = "seqpac",
                            parse="default_neb")
      }
      if(trim=="default_illumina"){
        counts<-make_counts(input=outpath,
                            trimming="seqpac",
                            parse="default_illumina")
      }
    }
    
  }
  if(is.null(lanes)){
    
    if(!is.null(trim)){
      cat("Will now merge lanes...")
      if(trim=="default_neb"){
        counts<-make_counts(input=inpath,
                            trimming="seqpac",
                            parse="default_neb")
      }
      if(trim=="default_illumina"){
        counts<-make_counts(input=inpath,
                            trimming="seqpac",
                            parse="default_illumina")
      }
    }
    
    if(is.null(trim)){
      counts<-make_counts(input=inpath)
    }
    
  }
  pheno <- data.frame(row.names = colnames(counts$counts),
                      Sample_ID= colnames(counts$counts))
  pheno<-make_pheno(pheno=pheno, progress_report=counts$progress_report,
                    counts=counts$counts)
  pac <- make_PAC(counts, pheno)
  
  closeAllConnections()
  return(pac)
}
