#' Advanced sequence mapping of a PAC object 
#'
#' \code{PAC_mapper} Map sequences against a small reference.
#'
#' Given a PAC object and the path to a fasta reference file, this function will
#' map sequences in the PAC using a 'backdoor' into the reanno workflow.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object.
#'
#' @param input Character indicating the path to the fasta (.fa) reference file or
#'   a DNAStringSet with already loaded reference sequences. If a Bowtie index
#'   is missing for the reference, PAC_mapper will attempt to temporarily
#'   generate such index automatically. Thus, large references are
#'   discouraged. Instead, we suggest you use the original reanno workflow 
#'   for large references.
#'   
#' @param mismatches Integer indicating the number of mismatches that should be
#'   allowed in the mapping.
#'
#' @param N_up Character indicating a sequence that should be added to the
#'   reference at the 5' end prior to mapping. A wild card nucleotides "NNN"
#'   (any of C, T, G, A) can for example be added for mapping non-perfect
#'   reference hits. No nucleotides are added by default.
#'   
#' @param N_down Character. Same as N_up but indicating a sequence that should
#'   be added to the reference at the 3' end. Useful for tRNA analysis where the
#'   reference do not contain pre-processed tRNA. Setting N_down="NNN" or "CCA"
#'   (in many species CCA is added to mature tRNA) will allow mapping against
#'   the mature tRNA. No nucleotides are added by default.
#'   
#' @param threads Integer indicating the number of parallel processes that
#'   should be used.
#'
#' @param report_string Logical whether an alignment string that shows where 
#'   sequences align against the reference in a character format. Works well with
#'   tRNA, but makes the Alignments object difficult to work with when longer
#'   references are used (default=FALSE).
#'
#' @param multi Character indicating how to deal with multimapping. If
#'   \code{multi="keep"}, query sequences that maps multiple times to the same
#'   reference sequence will be reported >1 times in the output (indicated by
#'   .1, .2, .3 etc. in the reported sequence name). If \code{multi="remove"}
#'   (default), then all multimapping sequences will be removed, resulting in 1
#'   row for each query sequence that maps to the target reference sequence. The
#'   function will always give a warning if a query sequence maps to multiple
#'   sites within a reference sequence. However, this function discriminate
#'   multimapping only within a reference sequence. Thus, if the fasta input
#'   contains multiple reference sequences, a query sequence may be reported in
#'   multiple references sequences.
#'   
#' @param override Logical whether or not the map_reanno function should prompt
#'   you for a question if there are files in the temporary path. As default,
#'   override=FALSE will prevent deleting large files by accident, but requires
#'   an interactive R session. Setting override=TRUE may solve non-interactive
#'   problems.
#'   
#'   
#' @return Stacked list, where each object on the highest level contains:
#'                    (Object 1) Reference name and sequence. 
#'                    (Object 2) Data.frame showing the mapping results of
#'                               each query sequence that mapped to Object 1.
#'
#' @examples
#' 
#'###########################################################
#'### Simple example of how to use PAC_mapper 
#' # Note: More details, see vignette and manuals.)
#' # Also see: ?map_rangetype, ?tRNA_class or ?PAC_trna, ?PAC_covplot 
#' # for more examples on how to use PAC_mapper.
#' 
#' ## Load PAC-object, make summaries and extract rRNA and tRNA
#'  load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                    package = "seqpac", mustWork = TRUE))
#' 
#' pac <- PAC_summary(pac, norm = "cpm", type = "means", 
#'                    pheno_target=list("stage", unique(pheno(pac)$stage)))
#'                    
#' pac_rRNA <- PAC_filter(pac, anno_target = list("Biotypes_mis0", "rRNA"))
#'
#' ## Give paths to a fasta reference (with or without bowtie index)
#' #  (Here we use an rRNA/tRNA fasta included in seqpac) 
#' 
#' ref_rRNA <- system.file("extdata/rrna", "rRNA.fa", 
#'                          package = "seqpac", mustWork = TRUE)
#'                          
#'
#' ## Map using PAC-mapper
#' map_rRNA <- PAC_mapper(pac_rRNA, mismatches=0, 
#'                         threads=1, input=ref_rRNA, override=TRUE)
#'
#' @export

PAC_mapper <- function(PAC, input, mismatches=0, multi="remove", 
                       threads=1, N_up="", N_down="", report_string=FALSE, 
                       override=FALSE){
  
  
  ## Setup
  j <- NULL
  ## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  ## Setup reference  
  if(methods::is(input, "DNAStringSet")){
    cat("\nImporting reference from DNAStringSet ...")
    full <- input
  }else{
    if(file.exists(input)){
      cat("\nReading reference from fasta file ...")
      full <- Biostrings::readDNAStringSet(input)
    }else{
      if(is.character(input)){
        cat("\nTry to import reference from character vector ...")
        full <- Biostrings::DNAStringSet(input)
      }else{
        stop("\nUnrecognizable reference format.",
             "\nPlease check your reference input.")
      }
    }
  }
  nams_full <- names(full) # Names are lost in the next step
  if(nchar(paste0(N_up, N_down)) >0){
    full <- Biostrings::DNAStringSet(paste(N_up, full, N_down, sep=""))
    names(full) <- nams_full
  }
  
  ## Setup temp folder and convert to windows format
  outpath <-  file.path(tempdir(), "seqpac")
  ref_path <-  file.path(tempdir(), "input","reference.fa")
  
  dir.create(outpath, showWarnings=FALSE, recursive = TRUE)
  dir.create(dirname(ref_path), showWarnings=FALSE, recursive = TRUE)
  
  Biostrings::writeXStringSet(full, filepath=ref_path, format="fasta")
  
  ## Make bowtie index if not available
  # If file input check bowtie index; save results in check_file
  check_file <- FALSE
  if(is.character(input)){
    if(file.exists(input)){
      if(!length(list.files(dirname(input), pattern=".ebwt"))>=2){
        check_file <- TRUE
      }
    }
  }
  if(check_file == FALSE){
    ref_path <- input 
    cat("\nBowtie indexes found. Will try to use them...")
  }else{
    if(nchar(paste0(N_up, N_down)) >0|
       methods::is(full, "DNAString")|
       check_file == TRUE){
      cat("\nNo bowtie indexes.")
      cat("\nWill try to reindex references ...")
      Rbowtie::bowtie_build(ref_path, outdir=dirname(ref_path), 
                            prefix = "reference", force = TRUE)  
    }
  }
  ## Make reanno object  
  map_reanno(PAC, input=list(reference=ref_path), output=outpath, 
             type="internal", threads=threads, mismatches=mismatches,  
             import="genome", keep_temp=FALSE, override=override)
  map <- make_reanno(outpath, PAC=PAC, mis_fasta_check = TRUE, output="list")
  stopifnot(length(map$Full_anno$mis0) == 1)
  
  ## Reorganize reanno object to a PAC_mapper object
  align_lst <- lapply(map$Full_anno, function(x){
    x <- x[[1]][!is.na(x[[1]]$ref_hits),]
    splt_x <- strsplit(x[[4]], "(?<=;\\+\\||;-\\|)", perl=TRUE)
    names(splt_x) <- x$seq
    lst_align <- lapply(splt_x, function(y){
      y <- gsub("\\|$", "", y)
      temp <- do.call("rbind", strsplit(y, ";"))
      nams <- temp[,1]
      start_align <- as.numeric(gsub("start=", "", temp[,2]))
      strnd_align <- temp[,3]
      strnd_align <- ifelse(strnd_align=="+", "sense", "antisense")
      df <- data.frame(ref_name = nams,
                       ref_strand = "*",
                       align_start = start_align,
                       align_strand = strnd_align)
    })
    df_align <- do.call("rbind", lst_align)
  })
  for(i in seq.int(length(align_lst))){
    align_lst[[i]]$seqs <- gsub("\\.\\d+", "", rownames(align_lst[[i]]))
    align_lst[[i]]$mismatch <- names(align_lst)[i]
  }
  # Rbind and fix names  
  nam_mis <- paste(paste0(names(align_lst), "."), collapse="|")
  align <- do.call("rbind", align_lst)
  rownames(align) <-  gsub(nam_mis, "", rownames(align)) 
  align_splt <- split(align, align$ref_name)
  align <- lapply(align_splt, function(x){
    rownames(x) <- NULL
    dup_tab <- table(x$seqs)
    nam_multi <- names(dup_tab)[dup_tab > 1]
    # Sequences mapping multiple times are removed or kept
    if(length(nam_multi) > 0){
      if(multi=="remove"){
        warning("\nSome sequences mapped >1 to the same reference.",
                "\nSince multi='remove' these sequences will be removed:",
                immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        x <- x[!x$seqs %in% nam_multi,]
        rownames(x) <- x$seqs
      }
      if(multi=="keep"){
        warning("\nSome sequences mapped >1 to the same reference.",
                "\nSince multi='keep', these sequences will be represented",
                "\nmultiple times in the mapping (psuedoreplication):",
                immediate. = TRUE)
        print(x[x$seqs %in% nam_multi,])
        splt <- split(x, x$seqs)
        splt <- lapply(splt, function(y){
          rownames(y) <- NULL
          return(y)
        })
        x <- do.call("rbind", splt)
      }
    }else{
      rownames(x) <- x$seqs
    }
    ifelse(x$align_strand=="sense", "+", "-")
    df <- data.frame(Mismatch=gsub("mis", "", x$mismatch),
                     Strand= ifelse(x$align_strand=="sense", "+", "-"),
                     Align_start=x$align_start, 
                     Align_end=x$align_start+nchar(x$seqs)-1, 
                     Align_width=nchar(x$seqs))
    rownames(df) <- rownames(x)
    return(df)
  })
  
  # Fix bowtie names and match with original reference
  splt_nam <- strsplit(names(full), " ")
  splt_nam <- unlist(lapply(splt_nam, function(x){x[1]}))
  nam_match <- match(splt_nam, names(align))
  align_lst <- align[nam_match]
  stopifnot(identical(names(align_lst)[!is.na(names(align_lst))],  
                      splt_nam[!is.na(names(align_lst))]))
  
  # Add full length reference
  names(align_lst)[is.na(names(align_lst))] <- splt_nam[is.na(names(align_lst))]
  fin_lst <- list(NULL)
  for(i in seq.int(length(align_lst))){
    if(is.null(align_lst[[i]])){
      align_lst[[i]] <- data.frame(Mismatch="no_hits", Strand="no_hits", 
                                   Align_start="no_hits", 
                                   Align_end="no_hits", Align_width="no_hits")
    }
    fin_lst[[i]] <- list(Ref_seq=full[i], Alignments=align_lst[[i]])
    names(fin_lst)[i] <- names(align_lst)[i] 
  }
  
  # Generate alignment string
  doParallel::registerDoParallel(threads) 
  `%dopar%` <- foreach::`%dopar%`
  
  if(report_string==TRUE){
    if(multi=="keep"){
      warning("\nOption multi='keep' is not compatible with report_string=TRUE",
              "\nAlignment string will not be returned.")
    }else{
      ref_lgn <- lapply(fin_lst, function(x){Biostrings::width(x$Ref_seq)})
      ref_lgn  <- max(do.call("c", ref_lgn))
      if(ref_lgn>500){
        warning("\nOption report_string=TRUE is only compatible with",
                "\nreference < 500 nt. Alignment string will not be returned.")
      }else{
        fin_lst <- lapply(fin_lst, function(x){
          if(x$Alignments[1,1] =="no_hits"){
            x$Alignments <- cbind(x$Alignments, 
                                  data.frame(Align_string="no_hits"))
            x$Alignments <- data.frame(lapply(x$Alignments, as.character), 
                                       stringsAsFactors=FALSE)
            return(x)
          }else{
            ref <- x$Ref_seq
            algn <- x$Alignments 
            n_ref <- nchar(as.character(ref))
            algn_lst <- split(algn, factor(row.names(algn), 
                                           levels=row.names(algn)))
            positions_lst <- foreach::foreach(j=seq.int(length(algn_lst)), 
                                              .final=function(y){
                                                names(y) <- names(algn_lst)
                                                return(y)})  %dopar% {
                                                  ref <- ref
                                                  n_ref <- n_ref
                                                  sq <- rownames(algn_lst[[j]])
                                                  if(algn_lst[[j]]$Strand == "-"){
                                                    sq <- intToUtf8(rev(utf8ToInt(sq)))
                                                  }
                                                  algn_str <- paste(strrep("-", 
                                                                           times=(algn_lst[[j]]$Align_start)-1),sq,
                                                                    strrep("-", 
                                                                           times= n_ref-(algn_lst[[j]]$Align_end)),
                                                                    sep="")
                                                  return(algn_str)
                                                }
            df <- cbind(algn, 
                        data.frame(
                          Align_string=as.character(
                            paste(do.call("c", positions_lst))), 
                          stringsAsFactors=FALSE))
            return(list(Ref_seq=ref, Alignments=df))
          }
        })
      }
    }
  }
  doParallel::stopImplicitCluster()
  #Have not implemented the the map object as a class in other functions yet
  return(fin_lst)
}
