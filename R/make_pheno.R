#' Make Phenotype File
#'
#' \code{make_pheno} generates a Pheno data.frame object containing additional
#' metadata for unique samples.
#'
#' Given the path to a directory with a single csv file, the function will read
#' the file. Alternatively, a data.frame with sample information can be
#' provided. If provided with a count table, (see \code{\link{make_counts}}), it
#' will attempt to organize the row names according to the column names in
#' the counts table.
#'
#' @family PAC generation
#'
#' @seealso  \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'   
#' @param pheno Character vector with the path to a .csv file or a
#'   data.frame with a column named "Sample_ID".
#'
#' @param counts Data.frame object with the same column names as in
#'   Sample_ID column of the pheno.
#'   
#' @param progress_report Data.frame object with progress report.
#'
#' @return A Pheno data.frame compatible with \code{\link{make_PAC}}
#'
#' @examples
#' 
#' ### First make counts 
#' 
#' # Seqpac includes strongly down-sampled smallRNA fastq.
#' sys_path = system.file("extdata", package = "seqpac", mustWork = TRUE)
#' input <- list.files(path = sys_path, pattern = "fastq", all.files = FALSE,
#'                 full.names = TRUE)
#'
#' # Notice that make_counts will generate another temp folder, that will 
#' # be emptied on finalization. By setting save_temp=TRUE you may save the 
#' # content.  
#'  
#' counts  <- make_counts(input, threads=2, parse="default_neb",
#'                        trimming="seqpac", plot=TRUE,
#'                        evidence=c(experiment=2, sample=1))
#'
#' colnames(counts$counts)
#' 
#' ### Then generate a phenotype table with make_pheno
#'
#' #  Note:  'Sample_ID' column needs to be similar IDs as 
#' #          colnames in the counts table. You may also 
#' #          specify a path to a txt file.
#'
#' Sample_ID <- colnames(counts$counts)
#'
#' pheno <- data.frame(Sample_ID=Sample_ID,
#'                        Treatment=c(rep("heat", times=1), 
#'                                    rep("control", times=2)),
#'                        Batch=rep(c("1", "2", "3"), times=1)) 
#' 
#' pheno <- make_pheno(pheno=pheno, progress_report=counts$progress_report, 
#'                      counts=counts$counts)
#'
#'
#' pheno 
#'
#' # Note that progress report from make_counts is added if you specify it  
#'      
#' ### Lastly combine into PAC
#' 
#' pac <- make_PAC(pheno=pheno, counts=counts$counts)
#'
#'
#' pac
#' names(pac)
#' 
#' # Note: a simple annotation table is added automatically.
#' head(anno(pac))
#'
#' # Clean up temp
#'closeAllConnections()
#'fls_temp  <- list.files(tempdir(), recursive=TRUE, full.names = TRUE)
#'file.remove(fls_temp, showWarnings=FALSE)
#'
#' @export
make_pheno<- function(pheno, counts=NULL, progress_report=NULL){
  
  ### Read using pheno.txt or data.frame 
  if(is.data.frame(pheno)){
    header <- which(grepl("^Sample_ID|^sample_ID|^Sample_id|^sample_id", 
                          colnames(pheno)))
    if(!length(header) == 1){
      stop(
        "\nCannot find column named 'Sample_ID'",
        "\nor you have >1 columns named 'Sample_ID'")
    }
    colnames(pheno)[header] <- "Sample_ID"
    
  }else{
    lines <- readLines(pheno, n=20)
    header <- which(grepl("Sample_ID|sample_ID|Sample_id|sample_id", lines))
    if(!length(header) == 1){
      stop("\nCannot find comma seperated header with first column", 
           "\nnamed 'Sample_ID', or you have >1 columns named 'Sample_ID'")
    }
    # Check separation by comparing1st and 2nd row (should have the same)
    head_1 <- stringr::str_count (lines[header], ",")
    row_1 <- stringr::str_count (lines[header+1], ",")
    head_2 <- stringr::str_count (lines[header], ";")
    row_2 <- stringr::str_count (lines[header+1], ";")
    head_3 <- stringr::str_count (lines[header], "\t")
    row_3 <- stringr::str_count (lines[header+1], "\t")
    
    if(row_1-head_1==0 & !head_1==0){
      pheno <- utils::read.delim(pheno,  header=TRUE, sep=",")
    }
    if(row_2-head_2==0 & !head_2==0){
      pheno <- utils::read.delim(pheno,  header=TRUE, sep=";")
    }
    if(row_3-head_3==0 & !head_3==0){
      pheno <- utils::read.delim(pheno,  header=TRUE, sep="\t")
    }
  }
  #Fixes bug in read.delim that messing up 1st ID_column
  id_col <- grepl("Sample_ID|sample_ID|Sample_id|sample_id", colnames(pheno))
  if(sum(id_col) >1){
    warning("There where multiple Sample_ID columns, will use the first one.")
  }
  colnames(pheno)[id_col] <- paste0("Sample_ID",seq_along(id_col[id_col]))
  colnames(pheno) <- gsub("Sample_ID1","Sample_ID", colnames(pheno))
  pheno$Sample_ID <- as.character(gsub("-", "_", 
                                       as.character(pheno$Sample_ID)))
  
  ## Order as counts using grepl
  if(!is.null(counts)){
    # Search column names in pheno$Sample_ID
    typ <- "pheno"
    ord <- unlist(lapply(as.list(colnames(counts)), function(x){
      logi_colnam <- grepl(x, pheno$Sample_ID)
      if(sum(logi_colnam)>1){
        stop("Input pheno does not have unique Sample_IDs")
      }
      if(sum(logi_colnam) == 0){
        return(0)
      }else{
        return(which(logi_colnam))
      }
    }))
    # If n matches does not match n pheno do the other way around
    # Save what type you have done until later
    if(!length(ord[!ord==0])==nrow(pheno)){
      typ <- "counts"
      ord <- unlist(lapply(as.list(pheno$Sample_ID), function(x){
        logi_colnam <- grepl(x,  colnames(counts))
        if(sum(logi_colnam)>1){
          stop("Input pheno does not have unique Sample_IDs")
        }
        return(which(logi_colnam))
      }))
    }
    
    # Double check everything and report missing samples from typ=pheno 
    if(any(duplicated(ord[!ord==0]))){
      stop("Sample_IDs were not unique in pheno input.")
    }
    if(!length(ord[!ord==0])==nrow(pheno)){
      stop(
        "\nNot all Sample_ID in pheno were available in counts column names.",
        "\nDouble check your Sample_ID column in pheno input.") 
    }
    
    # Reorder pheno according to counts sample names
    # Don't forget what type
    df <- as.data.frame(matrix(NA, nrow=ncol(counts), ncol=ncol(pheno)))
    rownames(df) <- colnames(counts)
    colnames(df) <- colnames(pheno)
    
    if(typ=="pheno"){
      for(i in seq.int(nrow(df))){
        if(!ord[i] == 0){
          df[i,] <- as.character(t(pheno[ord[i],]))
        }
      }
    }
    if(typ=="counts"){
      for(i in seq.int(nrow(df))){
        df[ord[i],] <- as.character(t(pheno[i,]))
      }
    }
    pheno <- df
    
    # Report outcome
    stopifnot(identical(rownames(pheno), colnames(counts)))
    logi_miss <- ord %in% 0
    print_df <- data.frame(pheno=as.character(pheno$Sample_ID), 
                           counts=colnames(counts))
    if(any(logi_miss)){
      warning(" Not all samples in counts were represented in pheno input.",
              "\n These will have 'NA' in pheno.")
    }
    cat("Of", length(colnames(counts)), 
        "sample names in counts,", 
        sum(!logi_miss), 
        "were found in pheno input.\n")
    if(!sum(print_df$pheno %in% print_df$counts) == length(print_df$counts)){
      cat("\nNote, partial name matching was done.")
    }
    cat("\n")
    print(print_df)
  }else{
    warning("\nNo count table was specified. Final Pheno will be unordered!\n")
  }
  
  ## Add progress report
  if(!is.null(progress_report)){
    cat("\nProgress report was specified, will attempt to match rownames...\n")
    
    # Fix sample names (illumina automatically exchanges "-" for "_"
    if(!sum(rownames(progress_report) %in% rownames(pheno)) == nrow(pheno)){
      rownames(progress_report) <- gsub("-", "_", rownames(progress_report))
    }
    prog_ord <- unlist(lapply(as.list(rownames(progress_report)), function(x){
      logi_colnam <- grepl(x, rownames(pheno))
      if(sum(logi_colnam)>1){
        stop("Input pheno does not have unique Sample_IDs")
      }
      return(which(logi_colnam))
    }))
    if(any(duplicated(prog_ord))){
      stop("Sample_IDs were not unique in pheno input.")
    }
    progress_report <- progress_report[match(rownames(pheno), 
                                             rownames(progress_report)),]
    stopifnot(identical(rownames(pheno), rownames(progress_report)))
    pheno <- cbind(pheno, progress_report)
    cat("\n")
  }else{
    warning("\nNo progress report was specified.",
            "\nThis will not be present in final Pheno.\n")
  }
  cat("Done!\n")
  return(pheno)
}




