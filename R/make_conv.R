#' Create reference genome conversion table 
#'
#' Generates a table that can be used for converting seqnames across referens
#' genomes.
#' 
#' Given the file paths to 2-3 fasta reference genomes, this function will
#' quickly match the the sequences between these reference genomes and return
#' seqnames of each matching sequence. Thus, given the paths to the fasta for
#' the Ensembl, UCSC and NCBI assemblies of a given genome version (e.g. hg38 or
#' dm6), will result in a name conversion table between these databases. 
#' 
#' Only perfect matches or no matches will be reported. Thus, in case an entire
#' sequence is missing (e.g. a sex chromosomes), a table will be returned with a
#' warning, but if a sequence is partly missing (e.g. only half a chromosome),
#' then an error is returned. 
#' 
#' Sequence matching is done using md5 hashes, which dramatically increases the
#' speed for perfect matches. Matching will always be done from reference A
#' against the other references. Thus, if reference A (=1st reference in
#' reference_list) contains less sequences than the other references, only the
#' reference A sequences will be reported in the output, having the same order
#' as in reference A.
#' 
#' @seealso \url{https://github.com/OestLab/seqpac} for updates.
#'
#' @family PAC reannotation
#'
#' @param reference_list List containing 2-3 file paths (character strings) to
#'   references in fasta format. The names of the list items will be imported
#'   into the final conversion table. The function may also take each reference
#'   path as single character strings (ref_path_A, ref_path_B, ref_path_C), but
#'   then any user-defined reference name will be lost (see example below).
#'   
#' @param ref_path_A Character string with path to 1st fasta reference (leading)
#' @param ref_path_B Character string with path to 2nd fasta reference
#' @param ref_path_C Character string with path to 3rd fasta reference
#'      
#' @param output Character indicating what type of output to provide. As default
#'   output="tibble" will result in a table in the tibble format as described in
#'   tibble/tidyverse packages. Anything else will result in a data frame.
#'   
#' @param skip_after Character or named list with character strings. If a
#'   character string, only the part of the seqnames prior to that string will
#'   be returned. Thus, if skip_after=" ", seqnames will be trimmed after first
#'   white space. If a list of character strings, then each string will will be
#'   applied to each provided reference. Thus, if skip_after=list(ensembl=" ",
#'   ucsc="," , ncbi=";"), then seqnames will be trimmed at white space for the
#'   reference named "ensembl", and at "," for the reference named "ucsc", etc.
#'   If skip_after=NULL then the full seqnames are returned.
#'    
#' @return A name conversion table either as a tibble or a data frame (see
#'   output)
#'   
#' @examples
#' 
#' ## Only for testing:
#' fasta_path <- system.file("extdata/trna", "tRNA.fa", 
#'                           package = "seqpac", mustWork = TRUE) 
#' ref1 <- Biostrings::readDNAStringSet(fasta_path) 
#' ref1 <- ref1[1:295]
#' sqnames <- do.call("rbind",(strsplit(names(ref1), "\\d chr")))[,2]
#' names(ref1) <- do.call("rbind",(strsplit(sqnames, " \\(")))[,1]
#' logi_dup <- duplicated(do.call("rbind", strsplit(names(ref1),"\\:"))[,1])
#' ref1 <- ref1[!logi_dup]
#' ref2 <- ref1
#' names(ref2) <- paste0("chr", names(ref2))
#' # Save new reference in temporary folder 
#' if(grepl("windows", .Platform$OS.type)){
#'   tmpdr <- paste0(tempdir(), "\\seqpac")
#' }else{
#'   tmpdr <- paste0(tempdir(), "/seqpac")}
#' dir.create(tmpdr, showWarnings=FALSE) 
#' ref_path1 <- paste0(tmpdr, "/ref1.fa")
#' ref_path2 <- paste0(tmpdr, "/ref2.fa") 
#' Biostrings::writeXStringSet(ref1, filepath=ref_path1, format="fasta")
#' Biostrings::writeXStringSet(ref2, filepath=ref_path2, format="fasta")
#' ref_list <- list(ensembl=ref_path1, ucsc=ref_path2)
#' conv_table <- make_conv(reference_list=ref_list)
#' conv_table
#'               
#'
#' ## The principles:
#' #
#' #ref_path_A <- "/some/path/to/ensembl.fa"
#' #ref_path_B <- "/some/path/to/ucsc.fa"
#' #ref_path_C <- "/some/path/to/refseq.fa"
#' #
#' #ref_list <- list(ensembl=ref_path_A, UCSC=ref_path_B, NCBI=ref_path_C)
#' #
#' ## Best (user defined names):
#' #conv_table <- make_conv(reference_list=ref_list) 
#' #
#' ## But also (no names)
#' #conv_table <- make_conv(ref_path_A, ref_path_B, ref_path_C)
#' #conv_table <- make_conv(ref_path_A, ref_path_C)
#' #
#' #
#' ## Make short names (skip everything after white space)
#' #conv_table <- make_conv(reference_list=reference_list, skip_after=" ") 
#' 
#' @export
make_conv <- function(reference_list=NULL,
                             ref_path_A=NULL, ref_path_B=NULL, ref_path_C=NULL,
                             output="tibble", skip_after=NULL){
  # Setup 
  sav_width <- list(NULL)
  if(is.null(reference_list)){
   reference_list <- list(refA=ref_path_A, refB=ref_path_B, refC=ref_path_C)
  }
  if(!is.null(skip_after)){
    if(!length(reference_list)== length(skip_after)){
      stop("\nReference length and skipafter length don't match.",
           "\nPlease provide input in the correct format.")
      }
  }
  if(length(reference_list)==2){
    reference_list <- c(reference_list , list(refC=ref_path_C))
  }
  
  # Generate md5 hash ref 1
  cat("\nReading reference A ...")
  ref1 <- Biostrings::readDNAStringSet(reference_list[[1]])
  stopifnot(!any(duplicated(names(ref1)))) 
  ref1_md5_vec <- unlist(lapply(as.list(ref1), function(x){
    digest::digest(paste(x) , algo="md5")
  }))
  sav_width[[1]] <- Biostrings::width(ref1)
  names(sav_width)[1] <- "ref1"
  rm(ref1)
  
  # Generate md5 hash ref 2
  cat("\nReading reference B ...")
  ref2 <- Biostrings::readDNAStringSet(reference_list[[2]])
  ref2_md5_lst <- lapply(as.list(ref2), function(x){
    digest::digest(paste(x) , algo="md5")
  })
  sav_width[[2]] <- Biostrings::width(ref2)
  names(sav_width)[2] <- "ref2"
  rm(ref2)
  
   # Generate md5 hash ref 3
  if(!is.null(reference_list[[3]])){
      cat("\nReading reference C ...")
      ref3 <- Biostrings::readDNAStringSet(reference_list[[3]])
      ref3_md5_lst <- lapply(as.list(ref3), function(x){
      digest::digest(paste(x) , algo="md5")
      })
      sav_width[[3]] <- Biostrings::width(ref3)
      names(sav_width)[3] <- "ref3"
      rm(ref3)
  }
  
  ## Ref1 vs Ref2
  # Check md5 between ref1 and ref2
  cat("\nMatching sequences in reference A vs B ...")
  hits1 <- character()
  for(i in seq.int(length(ref1_md5_vec))){
    frst_logi <- unlist(lapply(ref2_md5_lst, function(y){
          identical(paste(ref1_md5_vec[i]), paste(y))
       }))
    hits1 <- c(hits1, paste(which(frst_logi), collapse=";"))
   }
 
  # Check hits for duplicates and missing hits
    if(any(duplicated(hits1[!nchar(hits1) == 0]))){
      stop("\nThere were two or more identical reference ",
           "\nsequences having different names!")
    }
    if(any(grepl(";", hits1))){
      stop("\nThere were two or more identical reference ",
           "\nsequences having different names!")
    }
    if(sum(c(length(hits1) == length(ref1_md5_vec),
             length(hits1) == length(ref2_md5_lst)))==0){
      warning("\nNeither reference obtained full alignment.",
              "\nOutput will not report a complete conversion.", 
              immediate. = FALSE)
    }
    if(sum(c(length(hits1) == length(ref1_md5_vec), 
             length(hits1) == length(ref2_md5_lst)))==1){
      warning("\nOne of the references was not fully aligned.",
              "\nOutput will not report a complete conversion.", 
              immediate. = FALSE)
    }
  
  
  ## Ref1 vs Ref3
  if(!is.null(reference_list[[3]])){
    # Check md5 between ref1 and ref3
    cat("\nMatching sequences in reference A vs C ...")
    hits2 <- character()
    for(i in seq.int(length(ref1_md5_vec))){
      frst_logi <- unlist(lapply(ref3_md5_lst, function(y){
            identical(paste(ref1_md5_vec[i]), paste(y))
         }))
      hits2 <- c(hits2, paste(which(frst_logi), collapse=";"))
     }
   
    # Check hits for duplicates and missing hits
      if(any(duplicated(hits2[!nchar(hits2) == 0]))){
        stop("\nThere were two or more identical reference ",
             "\nsequences having different names!")
      }
      if(any(grepl(";", hits2))){
        stop("\nThere were two or more identical reference ",
             "\nsequences having different names!")
      }
      if(sum(c(length(hits2) == length(ref1_md5_vec),
               length(hits2) == length(ref3_md5_lst)))==0){
        warning("\nNeither reference obtained full alignment.",
                "\nOutput will not report a complete conversion.", 
                immediate. = FALSE)
      }
      if(sum(c(length(hits2) == length(ref1_md5_vec), 
               length(hits2) == length(ref3_md5_lst)))==1){
        warning("\nOne of the references was not fully aligned.",
                "\nOutput will not report a complete conversion.", 
                immediate. = FALSE)
      }
  }
  
  
  
  # Generate conversion table 
  cat("\nGenerating final conversion table ...")
  df <- data.frame(matrix(NA, nrow=length(hits1), ncol=6))
  colnames(df) <- c("name_ref_A", "name_ref_B", "name_ref_C", 
                      "width_ref_A", "width_ref_B", "width_ref_C")
  for(i in seq.int(length(hits1))){
        df$name_ref_A[i] <- names(ref1_md5_vec)[i]
        df$width_ref_A[i] <- sav_width$ref1[i]
        if(nchar(hits1[i])==0){
          df$name_ref_B[i] <- "no_hit"
          df$width_ref_B[i] <- "no_hit"
        }else{
          df$name_ref_B[i] <- names(ref2_md5_lst)[as.numeric(hits1[i])]
          df$width_ref_B[i] <- sav_width$ref2[as.numeric(hits1[i])]
        }
        
        if(!is.null(reference_list[[3]])){
            if(nchar(hits2[i])==0){
              df$name_ref_C[i] <- "no_hit"
              df$width_ref_C[i] <- "no_hit"
            }else{
              df$name_ref_C[i] <- names(ref3_md5_lst)[as.numeric(hits2[i])]
              df$width_ref_C[i] <- sav_width$ref3[as.numeric(hits2[i])]
            }
        }
    }
  
  # Check widths
  wdths <- df[!df$width_ref_B == "no_hit",]
  logi1 <- !as.numeric(wdths$width_ref_A) - as.numeric(wdths$width_ref_B)==0
  if(any(logi1)){
    stop("\nDiffering sequence lengths between fasta references.",
         "\nPlease make sure that you use the same genome versions.")
  }
  
  if(!is.null(reference_list[[3]])){
      wdths <- df[!df$width_ref_C == "no_hit",]
      logi2 <- !as.numeric(wdths$width_ref_A) - as.numeric(wdths$width_ref_C)==0
      if(any(logi2)){
      stop("\nDiffering sequence lengths between fasta references.",
           "\nPlease make sure that you use the same genome versions.")
  }
  }
  df <- df[, !grepl("width", colnames(df)), drop=FALSE]
  
  # Fix column names
  colnames(df) <- gsub("ref_A", names(reference_list)[1], colnames(df))
  colnames(df) <- gsub("ref_B", names(reference_list)[2], colnames(df)) 
  colnames(df) <- gsub("ref_C", names(reference_list)[3], colnames(df)) 
   
  # Modify names by skip after
  if(!is.null(skip_after)){
    if(!is.list(skip_after)){
      if(length(skip_after)>1){
        skip_after  <- as.list(skip_after)
      }else{
        skip_after  <- as.list(rep(skip_after, times=length(reference_list)))
      }
     names(skip_after) <- names(reference_list)
    }
    
    nam_check <- grepl(paste0(names(skip_after), collapse="|"), colnames(df))
    if(!sum(nam_check) == length(skip_after)){
      stop("\nReference names and skipafter names don't match.",
           "\nPlease provide input in the correct format.")
    }
  
    for(i in seq.int(length(skip_after))){
      clmn <- which(grepl(names(skip_after)[i],colnames(df))) 
                     

        df[,clmn] <- do.call("rbind", strsplit(df[,clmn], 
                                             split=skip_after[i]))[,1]
    }
  }
  if(output=="tibble"){
    return(tibble::as_tibble(df))
  }else{
    return(df)
  }
  cat("\n")
}
