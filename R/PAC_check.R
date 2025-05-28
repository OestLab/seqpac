#' Checks PAC format
#'
#' \code{PAC_check} Quickly checks rownames in PAC.
#'
#' Given a PAC object the function will cross check the rownames and colnames
#' across Pheno, Anno and Counts.
#' 
#' @family PAC generation
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object.
#' 
#' @param sample_check Can be provided either as a data frame or the name of a
#'   data frame in PAC. Will check both the row and column names of the object
#'   for identical sample names in PAC. A simple character vector can also be
#'   provided.
#' 
#' @param seq_check Can be provided either as a data frame or the name of an
#'   data frame in PAC. Will check only the row names of the object for
#'   identical sequence names in PAC. A simple character vector can also be
#'   provided.
#' 
#' @return If PAC format is correct will return logical TRUE, else it will throw
#'   an error.
#'   
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' 
#' PAC_check(pac)  # TRUE
#' 
#' # Remove first row in Anno and check compatibility
#' new_Anno <- anno(pac)[-1,]
#' #PAC_check(pac, seq_check=new_Anno) # Error
#' 
#' # Add to pac an check
#' # anno(pac) <- new_Anno #error

#' 
#' @export

PAC_check <- function(PAC, sample_check=NULL, seq_check=NULL){
  
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  
  sampl <- paste0("Please make sure that all samples are represented, ",
                  "\nnamed and ordered the correct way in both tables.") 
  seqs <-  paste0("Please make sure that all sequences are represented, ",
                  "\nnamed and ordered the correct way in both tables.")
  if(!identical(rownames(PAC$Pheno), colnames(PAC$Counts))){
    stop("\nSample names in Pheno (row names) and Counts",
         "(column names) are not identical.\n", sampl)
  }
  if(!identical(rownames(PAC$Anno), rownames(PAC$Counts))){
    stop("\nSequence names in Anno (row names) and Counts",
         "(row names) are not identical.\n", seqs)
  }
  if(length(PAC$norm) > 0){ 
    norm_logi_row  <- any(unlist(lapply(PAC$norm, function(x){
      !identical(rownames(PAC$Counts), rownames(x))
    })))
    norm_logi_col  <- any(unlist(lapply(PAC$norm, function(x){
      !identical(colnames(PAC$Counts), colnames(x))
    })))
    if(norm_logi_row){
      stop("\nRow names in normalized tables of PAC$norm are not identical",
           "\nwith row names in Counts.")
    }
    if(norm_logi_col){
      stop("\nColumn names in normalized tables of PAC$norm are not identical",
           "\nwith column names in Counts.")
    }
  }
  if(length(PAC$summary) > 0){ 
    sum_logi_row  <- any(!unlist(
      lapply(lapply(PAC$summary, rownames), function(x){
        identical(rownames(PAC$Counts), x)
      })))
    if(sum_logi_row){
      stop("\nRow names in summary tables of PAC$summary are not identical",
           "\nwith row names in Counts.")
    }
  } 
  res <- NULL
  ### sample_check #############################################################
  if(!is.null(sample_check)){
    ## Provided table
    if(methods::is(sample_check, c("data.frame", "matrix"))){
      if(any(rownames(sample_check)[1] %in% rownames(PAC$Pheno))){
        if(!identical(rownames(PAC$Pheno), rownames(sample_check))){
          stop("\nSample (row) names in new table are not identical",
               "\nwith sample names in the main PAC.\n", sampl)
        }
      } 
      if(any(colnames(sample_check)[1] %in% rownames(PAC$Pheno))){
        if(!identical(rownames(PAC$Pheno), colnames(sample_check))){
          stop("\nSample (column) names in new table are not identical",
               "\nwith sample names in the main PAC.\n", sampl)
        }
      }
    }else{
      ## Provided name of table in PAC
      if(is.character(sample_check)){
        if(length(sample_check)==1){
          whr  <- which(unlist(lapply(PAC, function(x){
            if(!methods::is(x, "list")){return(FALSE)}
            if(methods::is(x, "list")){any(grepl(sample_check, names(x)))}})))
          if(length(whr)==1){
            if(sum(grepl(paste(rownames(PAC$Pheno)[seq.int(2)], collapse="|"), 
                         rownames(PAC[[whr]][[sample_check]])))==2){
              if(!identical(rownames(PAC$Pheno), 
                            rownames(PAC[[whr]][[sample_check]]))){
                stop(
                  "\nSample (row) names in new table are not identical",
                  "\nwith sample names in the main PAC.\n", sampl)
              }
            }else{
              if(sum(grepl(paste(rownames(PAC$Pheno)[seq.int(2)], collapse="|"), 
                           colnames(PAC[[whr]][[sample_check]])))==2){
                if(!identical(rownames(PAC$Pheno), 
                              colnames(PAC[[whr]][[sample_check]]))){
                  stop(
                    "\nSample (column) names in new table are not identical",
                    "\nwith sample names in the main PAC.\n", sampl)
                }
              }else{   
                stop(
                  "\nPAC_check did find a matching object name, but could not",
                  "\nfind the samples. If you meant to check sequences ",
                  "\nand not samples, provide a seq_check object instead.")
              }
            }
          }
          
        }else{
          
          # Provided names only       
          if(length(sample_check)>1){ 
            if(!identical(rownames(PAC$Pheno), sample_check)){
              stop(
                "\nSample names are not identical with sample names in",
                "\nthe main PAC.\n", sampl)
            }
          }
        }
        
      }else{
        stop(
          "\nPAC_check did not recognize the sample_check input.",
          "\nDid you provide a data.frame, character vector or the",
          "\nname of a PAC object?")
      }  
    }}
  ### seq_check #############################################################
  if(!is.null(seq_check)){
    ## Provided table
    if(methods::is(seq_check, c("data.frame", "matrix"))) {
      if(!identical(rownames(PAC$Anno), rownames(seq_check))){
        stop(
          "\nSequence (row) names in new table are not identical with",
          "\nsequence names in the main PAC.\n", sampl)
      }
    }else{
      ## Provided name of table in PAC
      if(is.character(seq_check)){
        if(length(seq_check)==1){
          whr  <- which(unlist(lapply(PAC, function(x){  
            if(!methods::is(x, "list")){return(FALSE)}
            if(methods::is(x, "list")){any(grepl(seq_check, names(x)))}})))
          if(length(whr)==1){
            if(!identical(rownames(PAC$Anno), 
                          rownames(PAC[[whr]][[seq_check]]))){
              stop(
                "\nSequence (row) names in new table are not identical",
                "\nwith sequence names in the main PAC.\n", sampl)
            }
          }
        }else{
          ## Provided names only
          if(length(seq_check)>1){ 
            if(!identical(rownames(PAC$Anno), seq_check)){
              stop(
                "\nSequence names are not identical with sequence names",
                "\nin the main PAC.\n", sampl)
            }
          }}
      }else{
        stop(
          "\nThe PAC_check did not recognize the sample_check input.",
          "\nDid you provide a data.frame, character vector or the",
          "\nname of a PAC object?")
      }  
    }}
  return(TRUE)
}
