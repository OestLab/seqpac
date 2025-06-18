#'Trims/filter fastq using external cutadapt/fastq_quality_filter
#'
#'\code{make_cuadapt} cutadapt/fastq_quality_filter
#'
#'Given a path to sequence files in fastq format this function will trim adaptor
#'and remove sequences with low quality.
#'
#'@family PAC generation
#'
#'@seealso \url{https://cutadapt.readthedocs.io/en/stable/} for download and
#'documentation on cutadapt.
#'\url{http://hannonlab.cshl.edu/fastx_toolkit/commandline.html} for download
#'and documentation on fastq_quality_filter. \url{https://github.com/OestLab/seqpac}
#'for updates on seqpac.
#'
#'@param input Character path to a directory containing input fastq-files. The
#'  script will recursively search this directory for the .fastq|.fastq.gz
#'  extension.
#'
#'@param output Character path to the output directory where trimmed fastq files
#'  will be stored and temporary files will be generated.
#'
#'@param parse List with two character string expressions. The first will be
#'  parsed to cutadapt while the other is be parsed to fastq_quality_filter. If
#'  any is NULL, then the function will not pass the command and the trimming or
#'  filtering will not be applied. Thus, if parse = list(cutadapt=NULL,
#'  fastq_quality_filter="-q 20 -p 80"), then only the quality filter will be
#'  applied.
#'
#'@param threads Integer stating the number of parallel jobs. Note, that
#'  reading multiple fastq files drains memory fast, using up to 10Gb per fastq
#'  file. To avoid crashing the system due to memory shortage, make sure that
#'  each thread on the machine have at least 10 Gb of memory available, unless
#'  your fastq files are very small. Use \code{parallel::detectcores()} to see
#'  available threads on the machine.
#'
#'@return Externally the function will generate trimmed and/or quality filtered
#'fastq files in the output folder. Internally, a list of logs that can be used
#'to generate a progress report is returned.
#' 
#' @examples
#'  
#' ############################################################      
#' ### Principle of trimming using the make_cutadapt function
#' ### (Important: Need external installations of cutadapt 
#' ###  and fastq_quality_filter to work)
#' #  
#' #   input = "/some/path/to/input/folder"
#' #   output =  "/some/path/to/output/folder"
#' # 
#' ## Parse for make_cutadapt is a list of 2 character string expressions.
#' ## The first is parsed to cutadapt and the other to fastq_quality_filter 
#' ## For parallel processes '-j 1' is recommended since seqpac will   
#' ## parallelize across samples and not within.
#' ## Run system2("cutadapt -h", stdout=TRUE) and 
#' ## system("fastq_quality_filter -h", stdout=TRUE) 
#' ## for more options.
#' #   
#' ## String to parse to cutadapt:
#' # cut_prs <- paste0("-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAT",
#' #                    " --discard-untrimmed --nextseq-trim=20",
#' #                    " -O 10 -m 7 -M 70")
#' # 
#' ## Add string to parse to fastq_quality_filter:
#' #  parse = list(
#' #            cutadapt=cut_prs,
#' #            fastq_quality_filter="-q 20 -p 80")
#' #               
#' #  logs  <-  make_cutadapt(input, output, threads=8, parse=parse)
#' 
#' #' # Clean up temp
#'closeAllConnections()
#'fls_temp  <- list.files(tempdir(), recursive=TRUE, full.names = TRUE)
#'file.remove(fls_temp, showWarnings=FALSE)   
#'  
#' @export

make_cutadapt <- function(input, output, parse=NULL, threads=1){
  
  i <- NULL
  
  ## Setup
  cat("\nRunning make_cutadapt (R may stop responding) ...")
  cat("\n--- cutadapt and fastq_quality_filter must be correctly installed.")
  if(sum(!dir.exists(input))== length(input)){
    fls <- input
  }else{
    fls <- list.files(input, pattern ="fastq.gz\\>|\\.fastq\\>", 
                      full.names=TRUE, recursive=TRUE)
  }
  
  # Make dir
  if(!dir.exists(output)){
    logi_create <- try(dir.create(output), silent=TRUE)
    if(!is(logi_create, "try-error")){
      if(any(!logi_create)){
        warning("Was unable to create ", output, 
              "\nProbable reason: Permission denied")  
      
      }
    }
  }
  
  # Make output names and check output folder
  nam <- gsub("\\.gz$", "", basename(fls))
  nam <- gsub("\\.fastq$|\\.fq$|fastq$", "", nam)
  nam <- gsub("\\.$|\\.tar$", "", nam)
  nam_trim <- paste0(nam, ".trim.fastq.gz")
  # Check for output folder
  out_file <- file.path(output, nam_trim)
  out_dir <- list.files(output, pattern=nam_trim, recursive = FALSE)
  if(length(out_dir)>0){
    stop(
      "\n  Output trimmed fastq file names are identical to existing.",
      "\n  files in output:  ", out_dir, 
      "\n  Please move or delete the file in the output folder:\n  ", 
      output)
  }
  
  # Make sure parse is one line
  parse <- lapply(parse, function(x){gsub("[\r\n]", "", x)})
  
  ## cutadapt and fastq_quality_filter
  doParallel::registerDoParallel(threads) #Do not use parallel::makeClusters!!!
  `%dopar%` <- foreach::`%dopar%`
  prog_report <- foreach::foreach(i=seq.int(length(fls)), 
                                  .export= c("fls", "parse", "out_file", "nam"),
                                  .packages=c("ShortRead"), 
                                  .final = function(x){
                                    names(x) <- nam; return(x)}) %dopar% {
    log_lst <- list(NULL, NULL)
    names(log_lst) <- c("cutadapt", "fastq_quality_filter")
    spl_nam <- nam_trim[i]
    temp_out <- gsub("trim.fastq.gz$", "temp.fastq", out_file[i])
    if(!is.null(parse[[1]])){ 
      log_lst[[1]] <- system2(paste0("cutadapt ", parse[[1]], 
                                    " -o ", temp_out, " ", fls[i]), 
                              stdout=TRUE)
    }
    if(!is.null(parse[[2]])){ 
      log_lst[[2]] <- system2(paste0("fastq_quality_filter ", 
                                    parse[[2]], " -v -i ", temp_out, 
                                    " -o ", out_file[i], " -z"), 
                              stdout=TRUE)
    }
    
    log_in_trim <- grepl("Total reads processed", log_lst[[1]])
    log_w_adpt <- grepl("Reads with adapters", log_lst[[1]])
    log_shrt <- grepl("Reads that were too short", log_lst[[1]])
    log_lgn <- grepl("Reads that were too long", log_lst[[1]])
    log_out_trim <- grepl("Reads written \\(passing filters\\)", log_lst[[1]])
    log_in_q <- grepl("Input", log_lst[[2]])
    log_rm_q <- grepl("discarded", log_lst[[2]])
    log_out_q <- grepl("Output", log_lst[[2]])
    df <- data.frame(input=gsub(" |Total reads processed:|,", "", 
                                log_lst[[1]][log_in_trim]),
                     with_adapt= gsub(" |Reads with adapters:|,", "", 
                                      log_lst[[1]][log_w_adpt]),
                     too_short= gsub(" |Reads that were too short:|,", "", 
                                     log_lst[[1]][log_shrt]),
                     too_long= gsub(" |Reads that were too long:|,", "", 
                                    log_lst[[1]][log_lgn]),
                     out_adapt= gsub(" |Reads written \\(passing filters\\):|,",
                                     "", log_lst[[1]][log_out_trim]),
                     in_quality= gsub(" reads.|Input: |,", "", 
                                      log_lst[[2]][log_in_q]),
                     removed_quality= gsub("discarded | low-quality reads.", "",
                                           log_lst[[2]][log_rm_q]),
                     out_quality= gsub("Output: | reads.", "", 
                                       log_lst[[2]][log_out_q]),
                     stringsAsFactors=FALSE
    )
    df[1,] <-  gsub("\\(", " \\(", as.character(unlist(df[1,])))
    df[1,] <-  gsub("  ", " ", as.character(unlist(df[1,])))
    
    if(!is.null(parse[[1]])){  
      file.remove(temp_out)
    }
    return(df)
  }
  doParallel::stopImplicitCluster()
  prog_report <- do.call("rbind", prog_report)
  return(prog_report)
  gc(reset=TRUE)
}
