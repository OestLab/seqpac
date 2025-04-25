#' Imports annotation from reannotation mapping
#'
#' This function imports bowtie output and summarizes the content.
#'
#' Given the path to alignment outputs from bowtie, \code{\link{import_reanno}}
#' will attempt to read these files into R and generate a list of unsorted
#' data.frames where each row summarizes all annotations for a given sequence.
#' It is called by \code{\link{map_reanno}} to generate summarized data from
#' bowtie output saved as .Rdata files. These files can then be re-imported into
#' R where information is organized using \code{\link{make_reanno}} and finally
#' added to the annotation table (PAC$Anno) of the original PAC-list object
#' using \code{\link{add_reanno}}.
#'
#' @section Important: Re-annotation must have been done using Bowtie default
#'   output settings, where output files specifically have been named '.txt'.
#'   The function will not work with other formats such as SAM/BAM formats.
#'
#' @section Important: The basenames of the bowtie output will also be used to
#'   annotate. Example: bowtie -v 3 -a -f
#'   '<piRBase.fasta>' '<master_anno.fasta>' piRNA.txt Will annotate sequences
#'   as piRNA since the txt-file output was named 'piRNA'
#'
#' @family PAC reannotation
#'
#' @seealso  \url{http://bowtie-bio.sourceforge.net/index.shtml} for information
#'   about Bowtie and for Rbowtie:
#'   \url{https://www.bioconductor.org/packages/release/bioc/html/Rbowtie.html}.
#'   \url{https://github.com/Danis102} for updates on the current package.
#'
#' @param bowtie_path Path to a directory where bowtie output files can be
#'   found.
#'
#' @param report Character vector indicating what to report "minimum" or "full"
#'   (default="minimum").
#'
#' @param coord Logical whether or not mapping coordinates should be reported
#'   when report="full".
#'
#' @param reduce Character indicating a reference name (without file type
#'   extension) that should be exempted from report="full" and instead will
#'   generate a minimum report.
#'
#' @param threads Integer stating the number of parallel jobs.
#'
#' @return List of data frames with additional information from reannotation
#'   files generated with bowtie. If \emph{report="minimum"}, the function will
#'   report each hit only as the number of mismatches for a given reference
#'   file. If \emph{report="full"} the full name reported in the fasta file used
#'   as reference in the bowtie reannotation will be reported. If a reference
#'   name is specified in \emph{reduce}, this reference is excerpted from
#'   \emph{report="full"} and is instead reported as \emph{report="minimum"}.
#'
#'   Caution: Large references with lots of redundancy (such as pirBase in some
#'   species) will cause massive character strings if \emph{report="full"} with
#'   no restrictions. Specifying such references in
#'   \emph{reduce=<reference_names>} will circumvent this problem.
#'
#' @examples
#' 
#' ######################################################### 
#' ##### Test import_reanno
#'
#' ### First, if you haven't already generated Bowtie indexes for the included
#' # fasta reference you need to do so. If you are having problem see the small
#' # RNA guide (vignette) for more info.
#' 
#' ## tRNA reference:
#' trna_file <- system.file("extdata/trna", "tRNA.fa", 
#'                          package = "seqpac", mustWork = TRUE)
#' trna_dir<- gsub("tRNA.fa", "", trna_file)
#' 
#' if(!sum(stringr::str_count(list.files(trna_dir), ".ebwt")) ==6){
#'   Rbowtie::bowtie_build(trna_file, 
#'                         outdir=trna_dir, 
#'                         prefix="tRNA", force=TRUE)
#' }
#' ##  Then load a PAC-object and remove previous mapping from anno:
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                  package = "seqpac", mustWork = TRUE))
#' 
#' ref_paths <- list(trna= trna_file)
#' 
#' ##  You may add an output path of your choice, but here we use a temp folder:
#' output <- paste0(tempdir(),"/seqpac/test")
#' 
#' ##  Then map the PAC-object against the fasta references. Warning: if you use
#' # your own data, you may want to use override=FALSE, to avoid deleting previous
#' # mapping by mistake. keep_temp=TRUE can be used to run import_reanno
#' # independently.
#' 
#' map_reanno(pac, ref_paths=ref_paths, output_path=output,
#'            type="internal", mismatches=2,  import="biotype", 
#'            threads=2, keep_temp=TRUE, override=TRUE)
#' 
#' reanno1 <- import_reanno(output, report="Biotype",  threads=1)
#' 
#' @export

import_reanno <- function(bowtie_path, threads=1, coord=FALSE, 
                          report="minimum", reduce=NULL){
  s <- NULL
  base <- ".out$"
  files <- list.files(bowtie_path, pattern = base, full.names=TRUE)

  opt_sci <- options(scipen = 999)
  
  ## Check bowtie format (8 columns; "IIIIIII" present in column 6; 
  ## column 4 is an integer)
  row1 <- lapply(as.list(files), function(f){
    x <- try(utils::read.delim(f, nrows=1, header=FALSE), silent = TRUE)
    if(inherits(x, "try-error"))
      return(data.frame(V1="No_hits"))
    else
      return(x)
  })
  form_logi <- lapply(row1, function(x){
    if(!ncol(x) == 8){
      return(FALSE)
    }else{
      form_logi[[x]] <- return(sum(c(grepl("IIIIIII", 
                                           as.character(x[,6])), 
                                     is.integer(x[,4]))) == 2)
    }
  })
  form_logi <-  do.call("c", form_logi)
  no_hit <- do.call("c", lapply(row1, function(x){
    as.character(x[1,1]) == "No_hits"
  }))
  
  ## Give some feedback
  cat("\n|--- Found", sum(form_logi), 
      "bowtie file(s) with hits and", sum(no_hit), 
      "without.")
  
  ## Entering import loop
  data.table::setDTthreads(threads)
  bowtie_out_lst <- list(NA)
  for (k in seq.int(length(files))){
    cat(paste0("\n  |--- Import and reorganize ", basename(files)[k]))
    nam <- gsub(paste0(base), "",  basename(files)[k])
    
    ## Handle no hits    
    if(no_hit[k]){
      bowtie_out_lst[[k]] <- tibble::tibble(.id="No_hits", mis_n=NA,
                                            mis_where=NA, ref_hits=NA)
      names(bowtie_out_lst)[k] <- nam 
    }
    
    ## Import selected bowtie data
    if(form_logi[k]){
      
      if(report=="minimum"|nam %in% reduce){
        reprt <- "min"
        bow_out <- data.table::fread(files[k], header=FALSE, select = c(1,8), 
                                     data.table=TRUE, showProgress=FALSE)
      }else{
        reprt <- "full"
        if(coord==TRUE){bow_out <- data.table::fread(files[k], header=FALSE, 
                                                     select = c(3,4,2,1,8), 
                                                     data.table=TRUE, 
                                                     showProgress=FALSE)}
        if(coord==FALSE){bow_out <- data.table::fread(files[k], header=FALSE, 
                                                      select = c(3,2,1,8), 
                                                      data.table=TRUE, 
                                                      showProgress=FALSE)}
      }
      
      bow_out$V8 <- as.character(bow_out$V8)
      uni  <- unique(bow_out$V1)
      mis <- unique(bow_out$V8)
      n_mis <- stringr::str_count(mis, ":")
      n_mis <- as.character(unique(n_mis))
      if(is.na(n_mis)){n_mis <- 0}
      stopifnot(length(n_mis)=="1")
      n_mis <- paste0("mis", n_mis)
      
      ## Generate mini report from imported bowtie files
      if(reprt=="min"){
        cat(paste0("\n    |---> Generating minimum report ..."))
        bowtie_out_lst[[k]] <- tibble::tibble(.id=uni, mis_n=n_mis, 
                                              mis_where="mini_report", 
                                              ref_hits=nam)
        names(bowtie_out_lst)[k] <- nam
      }
      ## Compile full report with multithreading 
      if(reprt=="full"){
        cat("\n    |---> Generating full report (please wait)...")
        bow_splt <- split(bow_out, bow_out$V1)
        rm(bow_out)
        gc(reset=TRUE)
        
        # foreach combine every 100 instances
        chk_size <- ceiling(length(bow_splt)/100) 
        chnks1 <-as.integer(seq(from=1, to=length(bow_splt), by=chk_size))
        chnks2 <-as.integer(seq(from=0, to=length(bow_splt), by=chk_size))
        chnks2 <- c(chnks2[-1], length(bow_splt))
        chnks_rng <- list(chnks1, chnks2)
        
        # Do not use parallel::makeClusters!!!
        doParallel::registerDoParallel(threads) 
        `%dopar%` <- foreach::`%dopar%`
        bowtie_out_lst[[k]] <- foreach::foreach(s=seq.int(length(chnks_rng[[1]])), 
                .inorder = FALSE, 
                .combine = "rbind") %dopar% {
                  compile_lst <- lapply(
                    bow_splt[chnks_rng[[1]][s]:chnks_rng[[2]][s]], function(x){
                      # Fix neg strand mismatch   
                      new <- gsub("A", "t", x$V8[x$V2 == "-"])
                      new <- gsub("C", "g", new)
                      new <- gsub("G", "c", new)
                      new <- gsub("T", "a", new)
                      x$V8[x$V2 == "-"] <- toupper(new)
                      # Only report where mismatches occurs uniquely
                      uni_mis <- unique(x$V8)
                      uni_mis <- unique(
                        do.call("c", 
                                stringr::str_split(uni_mis, ",")))
                      uni_mis <- uni_mis[order(as.integer(gsub( ":.*$", "", 
                                                                uni_mis )))]
                      if(any(is.na(uni_mis))){uni_mis <- "mis0"}
                      uni_mis <- paste(uni_mis, collapse="|")
                      if(coord==TRUE){
                        x$V4 <- x$V4+1 # Fix bowtie coordinate shift
                        hits <- paste(unique(paste(x$V3, 
                                                   paste0("start=", 
                                                          x$V4), 
                                                   x$V2, sep=";")), 
                                      collapse="|")}
                      if(coord==FALSE){
                        strnd <- ifelse(x$V2=="+", "sense", "antisense")
                        hits <- paste(unique(paste(x$V3, strnd, sep=":")), 
                                      collapse="|")}
                      fin <- data.table::data.table(mis_n=n_mis, 
                                                    mis_where=uni_mis, 
                                                    ref_hits=hits)
                      return(fin)
                    })
                  bow_fin <- tibble::as_tibble(
                    data.table::rbindlist(compile_lst,idcol=TRUE))
                  return(bow_fin)
                }
        doParallel::stopImplicitCluster()
        names(bowtie_out_lst)[k] <- nam
      }
    }
    cat(paste0("\n    |---> ", nam, " done"))
  }
  
  on.exit(options(opt_sci)) 
  return(bowtie_out_lst)
}

