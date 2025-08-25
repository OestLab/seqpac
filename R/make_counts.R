#' Make a count table
#'
#'\code{make_counts} uses parallel processing to generate a count table.
#'
#' Given a path to fastq files this function performs low-level evidence filtering,
#' generating a counts table of sequences passing the filter and plots summary
#' statistics.
#' 
#' @family PAC generation
#'
#' @seealso  \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param input A path to a directory containing input fastq-files. The script
#'   will recursively search this directory for the .fastq|.fastq.gz extension.
#'
#' @param trimming Character indicating what type of trimming tool that should
#'   be used. If \code{trimming="seqpac"}, fastq files will be sent to the
#'   \code{\link{make_trim}} function in seqpac, while if
#'   \code{trimming="cutadapt"} fastq files will be sent to the
#'   \code{\link{make_cutadapt}} function. Note that \code{trimming="seqpac"}
#'   runs independently of external software, while \code{trimming="cutadapt"}
#'   is dependent on an externally installed version of cutadapt and
#'   fastq_quality_filter. Trimmed fastq files are stored temporarily in the
#'   systems default temporary folder. Please, run \code{\link{make_trim}} and
#'   \code{\link{make_cutadapt}} separately for permanent storage options, or
#'   set \code{save_temp=TRUE} to avoid that \code{make_counts} will delete all
#'   temporary files.  As default trimming=NULL, which indicates that input
#'   fastq files has already been trimmed.
#'     
#' @param parse Character strings defining the command that should be parsed to
#'   \code{\link{make_trim}} or \code{\link{make_cutadapt}}. This will allow
#'   you to customize your trimming according to 3' adapter sequence and
#'   platform standards etc. Please see examples below and the manuals for
#'   \code{\link{make_trim}} and \code{\link{make_cutadapt}} for more details.
#'   For convenience, \code{parse} also have two default mode for sRNA trimming,
#'   using Illumina and New England Biotype (neb) type small RNA adapters.
#'   \code{make_counts} will automatically print the exact setting for each
#'   default mode. Briefly, both modes involves polyG (NextSeq/NovaSeq) trimming
#'   and 3' adapter trimming, with a 0.1 tolerance for mismatch/indels. If
#'   parse="default_illumina", then the "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" 3'
#'   adapter is trimmed and untrimmed sequences are removed. If
#'   parse="default_neb", then "AGATCGGAAGAGCACACGTCTGAACTCCA" is trimmed and
#'   untrimmed sequences are removed. Removing untrimmed sequences is
#'   recommended for sRNA sequencing.
#'   
#' @param evidence Character vector with two inputs named 'experiment' and
#'   'sample' that controls the low-level evidence filter. Users may already at
#'   this point reduce the level of noise in the counts table by
#'   specifying the number of independent evidence that a specific sequence must
#'   have to be included. As default,
#'   \code{evidence=c(experiment=2, sample=1)} will include all sequences that
#'   have >=1 count in at least 2 independent fastq files. Thus 'experiment'
#'   controls the number of independent fastq evidence needed across the whole
#'   experiment. Note, however, that 'sample' does not control the number of
#'   counts needed in each sample. The evidence filter will always use >=1 count
#'   in X number of fastq files. Instead 'sample' controls when a sequence
#'   should be included despite not reaching the 'experiment' threshold. Thus if
#'   \code{evidence=c(experiment=2, sample=10)}, sequences that reach 10 counts
#'   in a single sample will also be included. If evidence=NULL all unique
#'   sequences will be saved. See 'examples' below for more examples.
#'   
#' @param optimize List that controls the behavior of make_counts,
#'   possibly enabling analysis on low-end computers on the expense of
#'   processing time. Has four items named: 
#'   list(on_disk=, chunk_size=, rm_wild=, seq_min=).
#'   1. on_disk: Logical whether or not some processes should be handled on-disk
#'   instead of in-memory. Important, this option will generate large temporary
#'   files. For a max fastq size of 1 GB (.gz compressed), keep at least 20 GB
#'   of free disk space. (default=FALSE).
#'   2. chunk_size: Integer, if set, determines whether reading and processing
#'   should be handled in chunks and how big those should be. If NULL, samples 
#'   are not handled in chunks. If handling a challenging dataset on low-end
#'   computers, a chunk_size of 50000000 is a good starting point.
#'   (default=NULL).
#'   3. rm_wild. Logical (panic filter). If evidence filter fails while
#'   on_disk=TRUE then rm_wild=TRUE will remove all sequences containing
#'   wild-card characters (N), either across the experiment or for individual
#'   samples. (default=FALSE).
#'   4. seq_min. Integer (panic filter).If evidence filter fails while
#'   on_disk=TRUE then seq_min will remove sequences shorter than seq_min,
#'   either across the experiment or for individual samples.(default=NULL)
#'   
#'   Note, Panic filters (2. and 3.) are applied at two places in the script as
#'   last resort if on_disk processing fails. First this may happen when unique
#'   sequences (UNI) fulfilling the evidence filter are extracted across the
#'   whole experiment. Here panic filters will be applied to all samples. The
#'   second place is when read sequences fulfilling the evidence filter are
#'   extracted from the trimmed fastq. Thus, the panic filter is applied to
#'   individual samples. To alert users, warnings are given in addition to notes
#'   in the panic_type column of the progress report. Here, "none" means no
#'   panic, "filt" means that panic filters were applied, and "chunk" means that
#'   loading of filtered fastq was forced into chunks. The last does not involve
#'   panic filters but helps loading very large fastqs on the expense of loading
#'   time.
#'   
#' @param plot Logical whether evidence plots should be printed and saved in the
#'   returned list (default=TRUE).
#'   
#' @param threads Integer stating the number of parallel jobs. Note, that
#'   reading multiple fastq files drains memory fast, using up to 10Gb per fastq
#'   file. To avoid crashing the system due to memory shortage, make sure that
#'   each thread on the machine have at least 10 Gb of memory available, unless
#'   your fastq files are very small. Use \code{parallel::detectcores()} to see
#'   available threads on the machine.
#'   
#' @param save_temp Logical whether temporary files (including trimmed fastq
#'   files) should be saved or not. Note, the function will print the path to
#'   the temprorary folder in the console.
#'   
#'   
#' @return 
#' A list containing three objects:
#' 
#'   1. counts (data frame) = count table. 
#'   
#'   2. progress_report = progress report from trimming and evidence filter.
#'   
#'   3. evidence_plots = bar graphs showing the impact of evidence filter. 
#'  
#' @examples
#' 
#'   
#' ############################################################ 
#' ### Seqpac fastq trimming with the make_counts function 
#' ### using default settings for NEBNext small RNA adaptor 
#' 
#' # Seqpac includes strongly down-sampled smallRNA fastq.
#' sys_path = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' input <- list.files(path = sys_path, pattern = "fastq", all.files = FALSE,
#'                 full.names = TRUE)
#'
#' # Now we can run make_counts
#' # Notice that make_counts will generate another temp folder, that will 
#' # be emptied on finalization. By setting save_temp=TRUE you may save the 
#' # content. You may also use your own trimming settings by using the parse 
#' # option. See ?make_trim for more details. 
#'  
#' counts  <- make_counts(input, threads=1, parse="default_neb",
#'                        trimming="seqpac", plot=TRUE,
#'                        evidence=c(experiment=2, sample=1))
#'
#' head(counts$counts)
#' head(counts$progress_report)
#' 
#' # Notice that that there are fewer unique sequences than number of reads 
#' # passing the filter. In normal, large, fastq we keep 95-99% of all reads 
#' # while filtering out 30-40% of the sequences that only appear in 1 fastq 
#' # Thus, the evidence filter may gain you performance in later analysis by 
#' # avoiding nonsense sequences. 
#'      
#'  
#'  #############################################################
#'  ### Lets change the evidence filter
#'  ###
#'  
#'  # 2 evidence over two indepenent samples, saving single sample 
#'  # sequences reaching 3 counts:
#'  
#'  test <- make_counts(input=input, trimming="seqpac", 
#'                      parse="default_neb",  
#'                      evidence=c(experiment=2, sample=3))
#'  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#'  test$counts[extras==1,]  # 6 single sample sequences reached 3 counts
#'  
#'  # 2 evidence over two independent samples, saving single sample 
#'  # sequences reaching 2 counts
#'   
#'  test <- make_counts(input=input,  trimming="seqpac", 
#'                      parse="default_neb",  
#'                      evidence=c(experiment=2, sample=2))
#'  extras <- apply(test$counts, 1, function(x){sum(!x==0)})
#'  test$counts[extras==1,] # 120 single sample sequences reached 2 counts
#'  
#'
#'  
#' @importFrom methods is
#' @export

make_counts <- function(input, trimming=NULL, threads=1, save_temp=FALSE,
                        plot=TRUE, parse="default_illumina", 
                        evidence=c(experiment=2, sample=1),
                        optimize=list(on_disk=FALSE, chunk_size=NULL,
                                      rm_wild=FALSE, seq_min=NULL)){
  
  anno <- value <- variable <- i <- NULL
  `%>%`<- dplyr::`%>%`
  `%dopar%` <- foreach::`%dopar%`
  opt_sci <- options()
  options(scipen=999)
  
  # Release optimization options
  chunk_size <- optimize[["chunk_size"]]
  on_disk <- optimize[["on_disk"]]
  rm_wild <- optimize[["rm_wild"]]
  seq_min <- optimize[["seq_min"]]
  #############################
  ###### Input fastq #######
  
  
  cat("Started at ", paste0(Sys.time()), "\n")
  Sys.sleep(1)
  
  gc(reset=TRUE)
  
  
  ## Read file system
  count_files <- list.files(input, pattern ="fastq.gz\\>|\\.fastq\\>", 
                    full.names=TRUE, recursive=TRUE, include.dirs = TRUE)
  if(length(count_files) == 0){
    count_files <- input
  }
  if(any(!file.exists(count_files))){
    stop("Something is wrong with input file(s)/path!")
  }
  count_files <- count_files[!grepl("Undetermined_", count_files)]
  count_files_nams <- basename(count_files)
  if(is.null(trimming)){
    cat("\nInput type was set to already trimmed fastq.")
  }else{
    cat("\nInput type was set to untrimmed fastq.")
  }
  cat("\nThe following fastq files were found in the path:\n")
  print(count_files)
  
  ## Setup temporary folder
  output <- file.path(tempdir(), "seqpac")
  if(!dir.exists(output)){
    dir.create(output)
  }
  if(!is.null(trimming)){
    if(dir.exists(output)){
      out_fls  <- list.files(output)
     logi_rm <- try(file.remove(out_fls), silent=TRUE) 
    }
    if(any(!logi_rm)){
      warning("Not all fastq was removed from ", output, 
              "\nProbable reason: Permission denied")
    }
    
    
    
    
  }else{
    trimming <- "trimmed"
  }
  
  ##########################################################
  ###### Trimming using cutadapt #######################
  if(trimming=="cutadapt"){
    
    cat("\nSeqpac trimming using the make_cutadapt function was specified.")
    cat("\n--- Temporary trimmed fastq are stored in: ", paste0(output))
    cat("\n--- Unless moved, these will be removed when system restarts.")
    cat("\n--- Please run 'make_cutadapt' for more options.")
    if(grepl("default", parse[[1]][1])){
      if(parse[[1]][1]=="default_neb"){
        cut_prse <- paste0("-j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCA",
                           " --discard-untrimmed --nextseq-trim=20", 
                           " -O 10 -m 14 -M 70 -e 0.1")
        parse <- list(cutadapt= cut_prse,
                      fastq_quality_filter="-q 20 -p 80")
      }
      if(parse[[1]][1]=="default_illumina"){ 
        
        cut_prse <- paste0("-j 1 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
                           " --discard-untrimmed --nextseq-trim=20",
                           " -O 10 -m 14 -M 70 -e 0.1")
        parse <- list(cutadapt=cut_prse,
                      fastq_quality_filter="-q 20 -p 80")
      }
      cat("\n--- Default trimming:\n   ")
      cat("\n---------------------------------------------")
      cat("\ncutadapt command:\n")
      print(parse$cutadapt)
      cat("\n---------------------------------------------")
      cat("\nfastq_quality_filter command:\n")
      print(parse$fastq_quality_filter)
      cat("\n\n")
    }
    
    ## Generate missing make_cutadapt inputs
    trim_input <- list("cutadapt", "fastq_quality_filter")
    names(trim_input) <- unlist(trim_input)
    trim_input <- lapply(trim_input, function(x){
      if(x %in% names(parse)){
        return(parse[x])
      }else{
        return(NULL)
      }
    })
    prog_report <- make_cutadapt(input=input, output=output, 
                                 threads=threads, parse=trim_input)
  }
  
  
  ##########################################################
  ###### Trimming using make_trim() #################
  if(trimming=="seqpac"){
    # Setup temporary file
    cat("\nSeqpac trimming using the make_trim function was specified.")
    cat("\n--- Temporary trimmed fastq are stored in: ", paste0(output))
    cat("\n--- Unless moved, these will be removed when system restarts.")
    cat("\n--- Please run 'make_trim' separately for more storage options.\n")
    if(grepl("default", parse[[1]][1])){
      if(parse[[1]][1]=="default_neb"){
        parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                      adapt_3_set=c(type="soft_rm", min=10, mismatch=0.1),
                      adapt_3="AGATCGGAAGAGCACACGTCTGAACTCCA",
                      seq_range=c(min=14, max=70), chunk_size=chunk_size,
                      quality=c(threshold=20, percent=0.8),
                      indels=TRUE, concat=12, check_mem=FALSE)
      }
      if(parse[[1]][1]=="default_illumina"){
        parse <- list(polyG=c(type="hard_trim", min=20, mismatch=0.1),
                     adapt_3_set=c(type="soft_rm", min=10, mismatch=0.1), 
                     adapt_3="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
                     seq_range=c(min=14, max=70), chunk_size=chunk_size,
                     quality=c(threshold=20, percent=0.8),
                     indels=TRUE, concat=12, check_mem=FALSE)
      }
      
      cat("\nDefault trimming:   \n")
      cat("----------------------------------------------")
      cat("\nadapt_3: \n")
      cat(parse$adapt_3, "\n")
      print(parse$adapt_3_set)
      cat("----------------------------------------------")
      cat("\npoly G: \n")
      print(parse$polyG)
      cat("----------------------------------------------")
      cat("\nsize filter:\n")
      print(parse$seq_range)
      cat("----------------------------------------------")
      cat("\nquality filter:\n")
      print(parse$quality)
      cat("----------------------------------------------\n")
    }
    
    ## Generate missing make_trim inputs for both default and parse
    trim_input <- list("polyG", "adapt_3", "adapt_3_set", "adapt_5", 
                       "adapt_5_set", "seq_range", "quality", "indels", 
                       "concat", "check_mem", "chunk_size")
    names(trim_input) <- unlist(trim_input)
    trim_input <- lapply(trim_input, function(x){
      if(x %in% names(parse)){
        return(parse[x])
      }else{
        if(x == "indels"){
          return(TRUE)
        }
        if(x == "concat"){
          return(12)
        }
        if(x == "check_mem"){
          return(TRUE)
        }
        if(x == "chunk_size"){
          return(chunk_size)
        }
      }
    })
    prog_report <- make_trim(input=input, output=output, threads=threads,
                             adapt_3=trim_input$adapt_3[[1]],
                             adapt_3_set=trim_input$adapt_3_set[[1]],
                             polyG=trim_input$polyG[[1]],
                             nucleotide_range = trim_input$seq_range[[1]],
                             quality=trim_input$quality[[1]],
                             concat=trim_input$concat[[1]],
                             indels=trim_input$indels[[1]],
                             check_mem=trim_input$check_mem[[1]],
                             chunk_size=trim_input$chunk_size[[1]]
    )
  }
  
  #########################################################
  ###### Read trimmed files and save sample evidence ######
  cat("\n\nIdentifying unique sequences in trimmed fastq files ...")
  if(!is.null(chunk_size)){
    cat("\nA chunk size of", chunk_size, "was specified...")
  }
  if(trimming=="trimmed"){
    fls <- count_files
  }
  if(trimming=="seqpac"){
    fls <- list.files(output, pattern="trim.fastq.gz$", 
                      full.names=TRUE, recursive=FALSE)
    fls <- fls[!grepl("REMOVED", fls)]
  }
  if(trimming=="cutadapt"){
    fls <- list.files(output, pattern="trim.fastq.gz$|trim.fastq$", 
                      full.names=TRUE, recursive=FALSE)
  }
  if(on_disk){  
    # Create path where sample uni/extra files end up prior to appending them.
    append_path <- file.path(output, "on_disk")
    if(dir.exists(append_path)){
      bad_fls <- list.files(append_path, full.names=TRUE)
      if(length(bad_fls)>0){
        file.remove(bad_fls)
      }
    }
    if(!dir.exists(append_path)){
      dir.create(append_path, recursive = TRUE)
    }
    cat("\nOn disk was requested. Progress for file append, see:\n",
        output)
  }
  names(fls) <- basename(fls)
  seq_lst <- lapply(fls, function(i){
    # Initiate
    fl<-i
    #fl<-fls[i]
    extras<-uni<-tab<-fl_nam_report<-NULL
    ## Read with loop for chunks
    if(!is.null(chunk_size)){
      #Uncompress the whole file will save time and space on disk
      # Otherwise vroom will generate a lot of files
      ugz_fl <- gsub(".gz|.gzip|.zip", "",basename(fl))
      ugz_path <- file.path(output, paste0("_tempugz_", ugz_fl))
      Sys.sleep(1)
      gc()
      gunzp_con <- R.utils::gunzip(fl, destname=ugz_path, 
                                   overwrite=TRUE, remove=FALSE)
      # Predict n chunks and end loop condition
      # ShortRead::countLines faster than R.util::countLines
      fq_lng <- ShortRead::countLines(ugz_path)
      
      tot_reads <- fq_lng/4
      tot_chunks <- fq_lng/chunk_size
      if(tot_chunks<=0){
        chunk_size <- fq_lng
        tot_chunks <- fq_lng/chunk_size
      }
      n_full_chunks <- floor(tot_chunks)
      perc_last_chunk <- tot_chunks-n_full_chunks
      last_chunk <- chunk_size*perc_last_chunk
      
      # Predict chunks
      chunks_fq <- seq(from =0 , to = fq_lng, by = chunk_size)
      
      # Run chunk for loop with last chunk in mind.
      # vroom_line is fastest but encounter problems at the end of very large 
      # fastq and generate big temp file. Eventually runs out of size.
      # readr::read_delim: slower but better at the end.
      # Reset `VROOM_CONNECTION_SIZE` to >131072*10000, helps to some degree
      # bigreadr::big_fread2 solves it completely but is relatively slow
      for (k in seq_along(chunks_fq)){
        start_chunk  <- chunks_fq[k]
        fq_left <- fq_lng - start_chunk
        
        if(fq_left >= chunk_size){
          fqtxt <- bigreadr::big_fread2(ugz_path, sep=NULL, select = 1,
                                        nrows=chunk_size, header=FALSE,
                                        skip=start_chunk, progress =FALSE,
                                        nThread=threads)
          
        }
        if(fq_left < chunk_size){
          fqtxt <- bigreadr::big_fread2(ugz_path, sep=NULL, select = 1, 
                                        header=FALSE, skip=start_chunk, 
                                        progress =FALSE, nThread=threads)
          Sys.sleep(1)
          unlink(ugz_path)
        }
        
        # Check which line contains seqs in chunk and extract seq lines
        # (every 4th; start depends on chunk_size; only the first 1000 rows)
        # Best to identify ID and strand lines in fastq, seqs are inbetween
        if(nrow(fqtxt) > 1000){
          samp_fq <- fqtxt[seq.int(1000),1]
        }else{
          samp_fq <- fqtxt[seq.int(length(fqtxt)),1]
        }
        # Find correct sample id and thereby correct id row.
        # Best to use multiple criteria to identify sample name from 1st 4 rows,
        # then use the name to find more rows
        id_sampl<- as.list(unlist(samp_fq)[seq.int(4)])
        # check common characters in id line
        find_id <- unlist(lapply(id_sampl, function(y) {
          sum(grepl("^@", y), grepl(":|.", y), 
              grepl(" ", y))==3}))
        # if that didnt work look for two lines with identical nchar =
        # seq and phred score lines
        if(!sum(find_id)==1){
          find_id <- lapply(id_sampl, function(y) {nchar(y)})
          strand_ln <- find_id==1
          S_Q_lns <- unlist(lapply(find_id, function(y){sum(y==find_id)==2}))
          find_id <- S_Q_lns + strand_ln
          find_id <- find_id == 0
        }
        if(!sum(find_id)==1){
          warning("Cannot identify fastq standard lines in chunks. Will start",
                  "\nfrom line 1 expecting chunk_size being a muliple of 4.")
          find_id <- c(TRUE,FALSE,FALSE,FALSE)
        }
        sampl_nam<- strsplit(as.character(id_sampl[unlist(find_id)]),
                             "\\.|:")[[1]][1]
        id_serie <- which(grepl(sampl_nam, samp_fq))
        
        strand_serie <- which(grepl("^[<\\+-]$", samp_fq))
        nt_serie <- strand_serie - id_serie
        nt_row <- sum(nt_serie)/length(nt_serie)
        if(!nt_row- round(nt_row)==0){
          warning("\nMay have detected unconventional fastq lines in chunk:\n",
                  basename(fl), "-", k)
        }
        nt_row <- round(nt_row)
        if(!nt_row %in% c(1,2,3,4)){ 
          stop("\nDetected severe unconventional fastq format in chunk:\n",
               basename(fl), "-", k)
        }
        
        seq_lines <- seq(from =nt_row , to = chunk_size-nt_row, by = 4)
        
        fqtxt <- fqtxt[seq_lines,]
        if(!is.null(evidence)){
          if(evidence[2] > 1){
            extras <- c(extras, table(fqtxt))
            extras <- tapply(extras, names(extras), sum)
          }
        }
        fqtxt <- unique(unlist(fqtxt, use.names = FALSE))
        fqtxt <- fqtxt[!is.na(fqtxt)]
        uni <- unique(c(uni, fqtxt))
        
        # Report chunk progress by writing a file
        if(!is.null(fl_nam_report)){
          file.remove(fl_nam_report)
        }
        fl_nam_report <- paste0("fastq", basename(fl), "_at_",
                                "_chunk_", k, "_of_", 
                                round(tot_chunks, digits=2))
        fl_nam_report <-  file.path(output, fl_nam_report)
        file.create(fl_nam_report)
        
        rm(fqtxt)
        gc(reset=TRUE)
      }
      # After chunk loop, extract extras. 
      # Uni will automatically be appended in loop but not extra.
      if(!is.null(evidence)){
        if(evidence[2] > 1){
          extras <- names(extras)[extras >= evidence[2]]
        }
      }
    }
    ## Without chunks, read whole file and extract unique seqs
    ##  put is.null(evidence) separate to prevent problems with NULL. 
    if(is.null(chunk_size)){
      if(is.null(evidence)){
        uni <- paste0(
          ShortRead::sread(
            ShortRead::readFastq(fl, withIds=FALSE)))
        uni <- unique(uni)
      }
      if(!is.null(evidence)){
        # For sample evidence >1 we need to extract extra using table
        if(evidence[2]>1){
          tab <- paste0(
            ShortRead::sread(
              ShortRead::readFastq(fl, withIds=FALSE)))
          tab <- table(tab)
          extras <- names(tab)[tab >= evidence[2]]
          uni <- names(tab)
          rm(tab)
        }
        # Only for evidence[2]=1, extract unique using ShortRead filters
        # Don't need table.
        if(evidence[2]==1){
          uni <- paste0(
            ShortRead::sread(
              ShortRead::readFastq(
                fl, withIds=FALSE, 
                filter=ShortRead::occurrenceFilter(duplicates="head", max=1)
              )))
        }
      }
    }
    ## Done reading and extract file
    # Return list depending on evidence and append on_disk
    if(!on_disk){
      if(is.null(evidence)){
        return(list(uni=uni, logi_extra=NULL))
      }
      if(!is.null(evidence)){
        if(evidence[2]==1){
          return(list(uni=uni, logi_extra=NULL))
        }
        if(evidence[2]>1){
          return(list(uni=uni, logi_extra=uni %in% extras))
        }
      }
    }
    # Too save memory append on_disk
    # make sure that parallel processes don't occupy the file
    # append multiple files afterward most secure
    if(on_disk){
      # First save seqs
      fl_nam_uni <- paste0(gsub(".fastq.gz", "", 
                                basename(fl)),"_APPEND_UNI.txt")
      fl_nam_uni <-  file.path(append_path, fl_nam_uni)
      # important; don't quote seqs (all must be handled equal downstream)
      readr::write_delim(as.data.frame(uni), 
                         fl_nam_uni,
                         append=FALSE,
                         col_names=FALSE,
                         quote="none")
      # Extra only if sample evidence >1
      if(evidence[2]>1){
        fl_nam_extra <- gsub("_APPEND_UNI.txt", 
                             "_APPEND_EXTRA.txt", fl_nam_uni)
        readr::write_delim(as.data.frame(extras), 
                           fl_nam_extra,
                           append = FALSE,
                           col_names=FALSE,
                           quote="none")
        
        return(list(uni=length(uni), logi_extra=length(extras)))
        rm(extras)
      }
      return(list(uni=length(uni), logi_extra=NULL))
    }
    # Clean up .Vroom and fastq temp files that may cause  
    # a full disk
    # Vroom uses lazy reading, file remove is only possible when 
    # the whole file has been read (may take time).
    rm(uni)
    gc()
    temp_fls <- list.files(tempdir(), full.names = TRUE, recursive=TRUE)
    vroom_fls <- temp_fls[grepl("vroom-",temp_fls)]
    if(length(vroom_fls)>0){unlink(vroom_fls)}
    tempugz_fls <- temp_fls[grepl("_tempugz_", temp_fls)]
    Sys.sleep(1)
    if(length(tempugz_fls)>0){unlink(tempugz_fls)}
    
    
  })
  Sys.sleep(1)
  
  gc(reset=TRUE)
  
  ###################################################################
  ### Append UNI and EXTRAS files if on_disk ####
  # If on_disk, append all files generated in foreach loop
  # Files may contain multiple of same seqs; evidence is used later
  
  cat("\nCompiling unique sequences ...")
  if(!on_disk){
    seqs <- unlist(lapply(seq_lst, function(x){x$uni}), use.names=FALSE)
  }
  if(on_disk){
    # First merge all sample UNI files to one file
    append_fls_uni<- list.files(append_path, 
                                pattern="_APPEND_UNI.txt", full.names=TRUE)
    fl_all_uni <-  file.path(append_path, "ALL_UNI.txt")
    file.append(fl_all_uni, append_fls_uni)
    file.remove(append_fls_uni)
    # If sample evidence2 > 1, do the same for extras
    if(evidence[2]>1){
      append_fls_extras<- list.files(append_path, 
                                     pattern="_APPEND_EXTRA.txt", 
                                     full.names=TRUE)
      fl_all_extras <-  file.path(append_path, "ALL_EXTRAS.txt")
      file.append(fl_all_extras, append_fls_extras)
      file.remove(append_fls_extras)
    }
    
    # Then read  merged UNI/EXTRA and reduce UNI over chunk while loop
    n_all_uni <- ShortRead::countLines(fl_all_uni)
    cat("\nOn_disk UNI file with fastq uniques had", n_all_uni, "sequences.")
    
    if(!is.null(chunk_size)){
      # Larger chunks will slow the counting; optimal with 10000000 
      if(chunk_size > 10000000){
        best_chk_size <- 10000000
      } else{
        best_chk_size <- chunk_size
      }
      if(n_all_uni <= best_chk_size){
        best_chk_size <- floor(n_all_uni/2)
      }
      Sys.sleep(1)
      
      gc()
      
      # Round seq.int( While-loop to convert UNI to fastq over chunks
      n_left <- n_all_uni        # n_left is all UNI from start
      start_chunk <- 0           # Start from 0 and add  best_chunk_size
      current_chunk <- 1         # Which chunk is the current?
      fl_all_good <-  file.path(append_path, "ALL_seqs_good.txt")
      fl_UNI_fastq <- file.path(append_path, "ALL_UNI.fq")
      keep_convert <- TRUE
      cat("\nProceeding with on disk fastq conversion of UNI file.")
      cat("\n")
      while(keep_convert){
        cat("\rChunk:", current_chunk, ". Seqs left to convert:", n_left)
        utils::flush.console()
        
        seqs_UNI <- bigreadr::big_fread2(fl_all_uni, sep=NULL, select = 1,
                                         header=FALSE, progress =FALSE,
                                         skip = start_chunk,
                                         nrows = best_chk_size,
                                         nThread=threads, nb_parts=1)
        seqs_UNI <- Biostrings::BStringSet(as.character(seqs_UNI[,1]))
        seqs_UNI <- seqs_UNI[order(seqs_UNI)]
        names(seqs_UNI) <- seq_along(seqs_UNI)
        Biostrings::writeXStringSet(seqs_UNI, fl_UNI_fastq, append=TRUE,
                                    format="fastq", compress=FALSE)
        
        current_chunk <- current_chunk +1
        start_chunk <- start_chunk + best_chk_size
        n_left <- n_all_uni - start_chunk
        if(start_chunk >= n_all_uni){
          keep_convert <- FALSE
        }
      } # Convert while loop ends
      
      #If not save, remove uni
      rm(seqs_UNI)
      if(!save_temp){
        unlink(fl_all_uni)
      }
      
      gc()
      
      # Round 2 Apply evidence[1] on disk using ShortRead::filterFastq
      if(evidence[2]>1){
        cat("\nApplying evidence filter ... ")
      }
      fl_filt_fastq <-  file.path(append_path, "Filt_UNI.fq")
      if(!is.null(evidence[1])){
        # ShortRead::filterFastq with ev1 function works!
        # gives warning when counting nts, probably because N wildcard present
        # evidence is applied anyway
        ev1_func <- function(x) {
          tab <- table(paste0(ShortRead::sread(x)))
          nams <-  names(tab)[tab >= evidence[1]]
          return(x[paste0(ShortRead::sread(x)) %in% nams])
        }
        err_ms <- try(ShortRead::filterFastq(files=fl_UNI_fastq, 
                                             destinations=fl_filt_fastq, 
                                             filter=ev1_func, 
                                             compress=FALSE), silent=TRUE)
        filt_uni <- try(ShortRead::readFastq(fl_filt_fastq, withIds=FALSE),
                        silent=TRUE)
        err_ms <- any(c(is(filt_uni, "try-error"), is(err_ms, "try-error")))
        if(err_ms){
          panic <- TRUE
        }else{
          panic <- FALSE #Needed later if individual samples fail
        }
        if(panic){
          warning("\nApplying evidence filter or loading filtered file failed.",
                  "\nThis may indicate a corrupt fastq or that many unique", 
                  "\nsequences are present even after evidence filter has", 
                  "\nbeen applied, which may exhaust memory. Will apply",
                  "\npanic filters to save the day...", immediate.=TRUE)
          
          fl_in <- fl_UNI_fastq
          fl_out <- fl_filt_fastq
          
          # seq_min likely leads to most reduction: apply first to save space 
          if(!is.null(seq_min)){
            min_func <- function(x) {
              x[nchar(paste0(ShortRead::sread(x))) >=seq_min]
            }
            ShortRead::filterFastq(files=fl_in, 
                                   destinations=fl_out, 
                                   filter=min_func,
                                   compress=FALSE)
          }
          if(rm_wild){
            # Update input/output only if a panic already has been applied 
            if(!is.null(seq_min)){
              fl_in <- fl_filt_fastq
              fl_out <- file.path(append_path, "Filt_UNI2.fq")
            }
            N_func <- function(x) {
               x[!grepl("N", ShortRead::sread(x))]
            }
            ShortRead::filterFastq(files=fl_in, 
                                   destinations=fl_out, 
                                   filter=N_func,
                                   compress=FALSE)
          }
          
          #Re-apply evidence filter after panic filters
          fl_in <- fl_out
          fl_out <- file.path(append_path, "Filt_UNI3.fq")
          #First clean up large UNI files
          if(grepl("Filt_UNI.", fl_in)){
            unlink(fl_UNI_fastq)
          }
          if(grepl("Filt_UNI2.", fl_in)){
            unlink(fl_filt_fastq)
          }
          message("\nNow reapplying evidence filters after panic filters...")
          ShortRead::filterFastq(files=fl_in, destinations=fl_out,
                                 filter=ev1_func, compress=FALSE)
          filt_uni <- try(ShortRead::readFastq(fl_filt_fastq, withIds=FALSE),
                          silent=TRUE)
          if(is(filt_uni, "try-error")){
            stop("\nReading sequences after evidence/panic filters failed.",
                 "\nThis indicates insufficient RAM memory, a very",
                 "\nlarge experiment, too many unique sequences possibly as a",
                 "\nconsequence of degraded samples. Try limit number input", 
                 "\nfastq or change the evidence/panic filters.")
          }
        } # Panic ends
        
        # Round 4 Make good seqs and append extras to good seqs 
        
        Sys.sleep(5)
        gc(reset=TRUE)
        filt_uni <- paste0(unique(ShortRead::sread(filt_uni)))
        # Read saved extras if evidence 2 and add to good
        if(evidence[2]>1){
          Sys.sleep(1) # will help in lowering mem imprint from previous vroom
          extras <- vroom::vroom_lines(fl_all_extras, num_threads =threads, 
                                       progress=FALSE, altrep = FALSE)
          filt_uni <- unique(c(filt_uni, extras))
          rm(extras)
        }
        readr::write_delim(as.data.frame(filt_uni),
                           fl_all_good,
                           append = FALSE,
                           col_names=FALSE,
                           quote="none")
      
        rm(filt_uni)
        Sys.sleep(5)
        gc(reset=TRUE)
      
        # Finally read GOOD.
        seqs_keep <- vroom::vroom_lines(fl_all_good, num_threads=threads, 
                                        progress=FALSE, altrep=FALSE)
        seqs_keep <- unique(seqs_keep)
        seqs_keep <- seqs_keep[order(seqs_keep, decreasing=FALSE)]

      Sys.sleep(10)
      gc(reset=TRUE)
      }#Evidence[1] done
    }#With chunks done
    
    #Without chunks read UNI immediately  
    if(is.null(chunk_size)){
      seqs <- readr::read_delim(fl_all_uni, delim="__", num_threads=threads, 
                                lazy =FALSE, col_names = FALSE,
                                progress=FALSE, show_col_types = FALSE)
      seqs <- unlist(seqs, use.names=FALSE)
    }
    if(evidence[2]>1){
        all_extras <- vroom::vroom_lines(fl_all_extras, num_threads =threads, 
                                     progress=FALSE, altrep = FALSE)
      }
    
    ## Save n uni seqs stats for on_disk, both with and without chunks
    n_uniseq <- unlist(lapply(seq_lst, function(x){
      x$uni}))
    
  } # on_disk ends
  
  #########################################################
  ###### Apply evidence filter and give info to user ######
  # Not on disk, but also on_disk here without chunks 
  # Evidence has not been extracted for without chunks 
  if(!on_disk|on_disk & is.null(chunk_size)){
    #All seqs regardless of input should be sorted  
    seqs <- seqs[order(seqs, decreasing=FALSE)]
    
    if(is.null(evidence)){
      seqs_keep <- unique(seqs)
    }else{
      if(evidence[1]==2){
        seqs <- seqs[duplicated(seqs)]
        seqs_keep <- unique(seqs)
        rm(seqs)
      }
      if(evidence[1]>2){
        seqs_keep <- table(seqs)
        rm(seqs)
        seqs_keep <- names(seqs_keep)[seqs_keep>=evidence[1]]
      }
      if(evidence[1]<2){
        seqs_keep <- unique(seqs)
        rm(seqs)
      }
    }
    
    ## Save n uni seqs stats for only not on_disk
    # on_disk without chunks, will receive this above
    if(!on_disk){
      n_uniseq <- unlist(lapply(seq_lst, function(x){
        length(x$uni)}))
    }
  } #End !on_disk | on_disk & is.null(chunk_size)
  
  ## evidence filter NULL message
  if(is.null(evidence)){
    cat(paste0("\nYou specified evidence=NULL. All sequences will be saved...",
               "\n(will take additional time and resources)"))
  }
  
  ## Sample (on disk sample evidence has already been extracted)
  if(!is.null(evidence)){
    cat(paste0("\nMaking a count table with sequences appearing in at least ", 
               evidence[1]," independent samples ..."))
    if(evidence[2]>1){
      cat("\nAlso saving sequences with less evidence but >= ", 
          evidence[2], " counts in a single sample.")
      cat("\n(Please set evidence$sample = 1 if this was not intended.)")
      
      if(!on_disk){
        lst_extras <- lapply(seq_lst, function(x){
          extras <- x$uni[x$logi_extra]
          n_extras <- sum(!extras %in% seqs_keep)
          return(list(extras=extras, n_extras=n_extras))
        })
        all_extras <- unlist(lapply(lst_extras, function(x){
          x$extras}), use.names=FALSE)
        rm(lst_extras)
      }
      # Last part needed for on_disk & is.null(chunk_size) too 
      seqs_keep <- unique(c(seqs_keep, all_extras))
      rm(all_extras)
      }
  }
  
  ###############################
  ###### Make count table ####### 
  # Don't closeAllCon in foreach loop! This will lead to failure.
  # Before and after instead!
  rm(seq_lst)
  Sys.sleep(10)
  gc(reset=TRUE)
  ## Prepare
  # If chunk_size or on_disk, user indicate low-end computer.
  # To avoid problems at the end here, reduce threads
  if(!is.null(chunk_size)|on_disk){
    threads <- floor(threads/2)
  }
  if(!is.null(chunk_size)&on_disk|threads<1){
    threads <- 1
  }
  doParallel::registerDoParallel(threads)
  `%dopar%` <- foreach::`%dopar%`
  check_done_path <- file.path(output, "check_fls_done.txt")
  # panic needs to start both inside and outside loop
  # unless returned (not needed), it must restart outside loop for !on_disk
  # on_disk it needs to restart inside loop or it won't be found
  panic <- "none"
 
  ## Collect selected reads from fls using foreach loop
  reads_lst <- foreach::foreach(i=seq_along(fls),.final = function(x){
   names(x) <- basename(fls); return(x)}) %dopar% {
      panic <- "none"
      stat_mat <- as.data.frame(matrix(NA, nrow=1, ncol=4))
      colnames(stat_mat) <- c("tot_reads", "reads_pass_evidence", 
                              "uniseqs_pass_evidence", "panic_type")
      fl <- fls[[i]]
      filt_err<-filt_err2<-reads_keep <- NULL
      # If on disk, extract seqs prior to reading the fastq using filterFastq
      # This will save some ram but will temporally add lot to temp folder
      if(on_disk){
        fq_lng <- ShortRead::countLines(fl)
        stat_mat$tot_reads <- fq_lng/4
        
        #Modify chunk_size if too much 
        #Add if chunk_size is missing (may happen if only on_disk is provided)
        best_chk_size <- 10000000
        if(!is.null(chunk_size)){
          if(chunk_size < 10000000){
            best_chk_size <- chunk_size
          }
        }
        if(fq_lng <= best_chk_size){
           best_chk_size <- floor(fq_lng/2)
        }
        
        # Fitler/panic functions
        extract_fun <- function(x) {
          x[paste0(ShortRead::sread(x)) %in% seqs_keep]
        }
        pfilt_fun <- function(seqs_in, fl_in, fl_out, extract_fun,
                              seq_min, rm_wild, threads){
          # First save seqs_keep
          fl_seq_sav<- file.path(tempfile(),"_temp_save_seqs.txt")
          readr::write_delim(as.data.frame(seqs_in),
                             fl_seq_sav,
                             append = FALSE,
                             col_names=FALSE,
                             quote="none")
          
          if(!is.null(seq_min)){
            seqs_in <- seqs_in[nchar(seqs_in) >=seq_min]
          }
          if(rm_wild){
            seqs_in <- seqs_in[!grepl("N", seqs_in)]
          }
          gc()
          filt_err<- ShortRead::filterFastq(files=fl, 
                                                 destinations=fl_temp, 
                                                 filter=extract_fun, 
                                                 compress=FALSE)
          # Must re-load seqs_keep
          seqs_keep <- vroom::vroom_lines(fl_seq_sav, num_threads=threads,
                                          progress=FALSE, altrep=FALSE)
          return(filt_err)
        }
        chunk_fun <- function(fl, chunk_size){
          Sys.sleep(10)
          gc()
          reads_keep <- NULL
          fq <- 1
          sampler <- ShortRead::FastqStreamer(fl_temp, chunk_size)
          while(length(fq)) {
            fq <- ShortRead::yield(sampler)
            reads_keep <-c(reads_keep,
                               table(paste0(ShortRead::sread(fq))))
            reads_keep <-tapply(reads_keep, names(reads_keep), sum)
          }
          rm(fq)
          close(sampler)
          n_reads <- as.data.frame(reads_keep)
          return(n_reads)
        } # chunk_func ends
        
        # Start error handling  
        filt_out <- NULL
        tmp_fln_nam<- paste0(i, "_temp", basename(fl),".fq")
        fl_temp <-  file.path(output, tmp_fln_nam)
        gc()
        # Apply panic in rounds, round1=filt, round2=chunk, round3=remove
        # Try filter, if error, apply panic filter to this sample only
        filt_out <-try(ShortRead::filterFastq(files=fl, destinations=fl_temp, 
                                              filter=extract_fun, 
                                              compress=FALSE), silent=TRUE)
        #If filter failed then apply pfilt_func  
        if(is(filt_out,"try-error")){
          warning("Failed to filter \n", fl, 
                  "\non_disk. Will try to save the day by applying panic",
                  "\nfilters to this sample. Important, before proceeding", 
                  "\nwith these counts you must exlude or harmonize this",
                  "\nsample with the rest of the experiment by filtering the",
                  "\nother samples.")
          panic <- "filt"
          filt_out <- try(pfilt_fun(seqs_in=seqs_keep, fl_in=fl, fl_out=fl_temp,
                                    extract_fun=extract_fun, threads=threads,
                                    seq_min=seq_min, rm_wild=rm_wild), 
                          silent=TRUE)
        }
        # If filter worked, try to read filtered temp file
        if(!is(filt_out,"try-error")){
          Sys.sleep(10)
          gc()
          reads_keep <- try(ShortRead::readFastq(fl_temp, withIds=FALSE),
                            silent=TRUE)
        }
        #If reading failed, then apply chunk_fun 
        if(is(reads_keep,"try-error")){
          Sys.sleep(10)
          gc()
          warning("Failed to read ", fl, "",
                  "\non_disk after evidence filter. Will try to save the", 
                  "\nday by reading this sample in chunks, which will take",
                  " additional time.")
          if(panic=="filt"){
            panic <- "filt_chunks"
          }else{
            panic <- "chunks"
          }
          reads_keep <- try(chunk_fun(fl=fl_temp, chunk_size=best_chk_size), 
                            silent=TRUE)
        }
        #If all fail, remove sample
        if(is(reads_keep,"try-error")){
          # Temporarily set seqs to reads to avoid problems later
          warning("\nFailed to filter/read", fl, "again.",
                  "\nWill exlude this sample from counts")
          panic <- "remove"
          reads_keep <- seqs_keep
        }
        gc()

        ## Clean up large temporary file
        #First try to remove file
        if(any(file.exists(fl_temp))){
          Sys.sleep(5)
          unlink(fl_temp)
        }
        # If it doesnt work, append file path to file 
        # Try to remove all finished files in every loop
        write(fl_temp, file=check_done_path, append=TRUE)
        done_fls<- readLines(check_done_path)
        if(any(file.exists(done_fls))){
          unlink(done_fls)
        }
        
        # Make counts of kept reads: n_reads
        if(panic=="none"){
          reads_keep <- paste0(ShortRead::sread(reads_keep))
          n_reads <- as.data.frame(table(reads_keep))
          stat_mat$panic_type <- "none"
        }
        #Panic chunk is already n_reads table (data.frame)
        if(panic=="chunks"|panic=="filt_chunks"){
          reads_keep$temp <- rownames(reads_keep)
          colnames(reads_keep) <- c("Freq", "reads_keep")
          rownames(reads_keep) <- NULL
          n_reads <- as.data.frame(reads_keep[,c(2,1)])
        }
        if(panic=="chunks"){
          stat_mat$panic_type <- "chunks"
        }
        if(panic=="filt_chunks"){
          stat_mat$panic_type <- "filt_chunks"
        }
        if(panic=="filt"){
          stat_mat$panic_type <- "filt"
        }
        if(panic=="remove"){
          dt["N"] <- NA
          stat_mat$panic_type <- "removed"
          stat_mat$reads_pass_evidence <- NA
          stat_mat$uniseqs_pass_evidence <- NA
        }
      rm(reads_keep)
      } #on_disk ends
      
      if(!on_disk){
        reads <- ShortRead::readFastq(fl, withIds=FALSE)
        reads <- paste0(ShortRead::sread(reads))
        stat_mat$tot_reads <- length(reads)
        stat_mat$panic_type <- "none"
        if(is.null(evidence)){
          reads_keep <- reads
        }
        if(!is.null(evidence)){
          reads_keep <- reads[reads %in% seqs_keep]
        }
        n_reads <- as.data.frame(table(reads_keep))
        rm(reads,reads_keep)
      } #not on_disk ends
      
      # To all, add stat prior to match; match introduces NA rows (0 counts)
      seqs_pass <- try(nrow(n_reads))
      if(is(seqs_pass, "try-error")){
        seqs_pass <-  length(n_reads) 
      }
      stat_mat$uniseqs_pass_evidence <- seqs_pass
      stat_mat$reads_pass_evidence <- sum(n_reads$Freq)
      n_reads <- n_reads[match(seqs_keep, n_reads$reads_keep),,drop=FALSE]
      n_reads$Freq[is.na(n_reads$Freq)] <- 0
      dt <- tibble::tibble(N=n_reads$Freq)
      rm(n_reads)
      return(list(counts=dt, stat=stat_mat))
      gc()
    }# foreach loop ends
  doParallel::stopImplicitCluster()
  gc(reset=TRUE)
  
  ## Finalize 
  cat("\nFinalizing at ", paste0(Sys.time()), "\n")
  names(reads_lst) <- gsub("_merge|\\.fastq.gz$|fastq.gz$|\\.fastq$", 
                           "", names(reads_lst))
  names(reads_lst) <- gsub("temp_|\\.trim", "", names(reads_lst)) 
  names(reads_lst) <- gsub("-", "_", names(reads_lst))
  
  ## Clean up temporary trimmed files
  if(!trimming=="trimmed"){
    if(save_temp==FALSE){
      fn <- list.files(output, full.names=TRUE, recursive=FALSE)
      if(any(file.exists(fn))){
        file.remove(fn)
      }
    }
  }
  # restore vroom env settings and clean up after vroom
  if(on_disk){
    gc(reset=TRUE)
    temp_fls <- list.files(tempdir(), full.names = TRUE)
    vroom_fls <- temp_fls[grepl("vroom-",temp_fls)]
    if(any(file.exists(vroom_fls))){
      unlink(vroom_fls)
    }
  }
  # Build df
  ordCount_df <- as.data.frame(do.call("cbind", lapply(reads_lst, function(x){
    x$counts})), stringsAsFactors = FALSE)
  colnames(ordCount_df) <- names(reads_lst)
  rownames(ordCount_df) <- seqs_keep
  
  if(trimming=="trimmed"){
    prog_report <- data.frame(
      trimming=as.character(rep(
        "no_trimming_was_done", times=length(reads_lst))), 
      stringsAsFactors=FALSE)
  }
  
  ordCount_df <- ordCount_df[order(rownames(ordCount_df)),]
  
  ###################################################
  ###### Add evidence stats to progress report ######
  stat_dt <- cbind(data.frame(uni_seqs=n_uniseq, stringsAsFactors =FALSE), 
                    do.call("rbind", lapply(reads_lst, function(x){x$stat})))
  stat_dt$reads_not_pass_evidence <- stat_dt$tot_reads - 
     stat_dt$reads_pass_evidence
  stat_dt$uniseqs_not_pass_evidence <- stat_dt$uni_seqs - 
     stat_dt$uniseqs_pass_evidence
  rownames(stat_dt) <- colnames(ordCount_df)
  stat_dt$sample <- colnames(ordCount_df) 
  if(max(nchar(stat_dt$sample))>15){
    stat_dt$sample <- paste0("fastq_", seq.int(length(stat_dt$sample)))
  }
  df_long <- reshape2::melt(stat_dt, id.vars="sample")
  df_long$value <- as.numeric(df_long$value)
  df_long$sample <- as.factor(df_long$sample)
  main_title <- "Reads passed that evidence"
  
  if(is.null(evidence)){
    evidence <- c(experiment=1, sample=1)
  }
  
  p_func <- function(df, y_lab, main_title, sub_title){
    df$variable <- factor(df$variable)
    ord <- c(which(grepl("_not_", levels(df$variable))), 
             which(!grepl("_not_", levels(df$variable)))) 
    df$variable <- factor(df$variable, 
                          levels=levels(df$variable)[ord]) 
    ggplot2::ggplot(df, ggplot2::aes(x=sample, y=value, fill=variable)) +
      ggplot2::geom_col(width=1, size=0.2, color="black") + 
      ggplot2::geom_hline(yintercept=0, col="#707177", cex=0.6) +
      ggplot2::geom_text(ggplot2::aes(label=value), 
                         position = ggplot2::position_stack(vjust = 0.5), 
                         angle = 90, color="white", size=4)+
      ggplot2::scale_fill_manual(name = NULL, labels = c("not pass", "pass"), 
                                 values=rev(c("#085B84", "#B20018")))+
      ggplot2::ylab(y_lab)+
      ggplot2::labs(title=main_title, subtitle=)+
      ggplot2::theme(
        plot.caption =  ggplot2::element_text(size=10, face= "bold"),
        axis.title.y = ggplot2::element_text(size=12, face= "bold"),
        axis.title.x = ggplot2::element_text(size=12, face="bold"),
        axis.text = ggplot2::element_text(size=10),
        axis.text.x = ggplot2::element_text(angle=45, hjust=1),
        panel.background = ggplot2::element_blank())
  }
  p1 <- p_func(df_long[df_long$variable %in% c("reads_pass_evidence", 
                                               "reads_not_pass_evidence"),], 
               y_lab="number of reads", 
               main_title=paste0("Sequences present in at least ", 
                                 evidence[1], " sample(s)"),
               sub_title=if(evidence[2]>1){paste0("(unless >", 
                                                  evidence[2], 
                                                  "reads in a single sample)")
               }else{
                 NULL
               })
  p2 <- p_func(df_long[df_long$variable %in% c("uniseqs_pass_evidence", 
                                               "uniseqs_not_pass_evidence"),], 
               y_lab="number of unique sequences", 
               main_title=paste0("Sequences present in at least ", 
                                 evidence[1], " sample(s)"),
               sub_title=if(evidence[2]>1){paste0("(unless >", evidence[2], 
                                                  "reads in a single sample)")
               }else{
                 NULL}
  )
  plt_lst <- list(reads=p1, uniseqs=p2)
  if(plot==TRUE){
    print(cowplot::plot_grid(plotlist=plt_lst, ncol=1, nrow=2))
  }
  if(plot==FALSE){
    plt_lst <- "Evidence plot was omitted by user input"
  }
  prog_report <- cbind(prog_report, stat_dt)
  options(opt_sci)
  return(list(counts=ordCount_df, progress_report=prog_report, 
              evidence_plots=plt_lst))
}
