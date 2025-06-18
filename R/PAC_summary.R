#' Summarizes PAC objects 
#'
#' \code{PAC_summary} summarizes data stored in a PAC object.
#'
#' Given a PAC object this function summarize data in counts(PAC) or in the norm
#' 'folder' according to a grouping column in pheno(PAC).
#' 
#' @family PAC analysis
#' 
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#' package.
#'
#' @param PAC PAC object containing a Pheno data.frame with samples as row
#'   names and a Counts table with raw counts. Optionally, the PAC object may
#'   contain normalized counts tables, saved in the norm list ('folder'). Such
#'   normalized table can be generated using the \code{\link{PAC_norm}}
#'   function.
#'   
#' @param norm Character indicating what type of data to be used. If 'counts',
#'   the raw counts in Counts will be used (default). Given any other value, the
#'   function will search for the value as a name on a data.frame stored in the
#'   normalized list-folder.
#'
#' @param type Character indicating what type of summary to be applied to the
#'   data. The function currently supports:
#'   type="means"         # Group means
#'   type="sd"            # Group standard deviation
#'   type="se"            # Group standard error of the mean
#'   type="log2FC"        # Group log2 fold changes against other groups
#'   type="log2FCgrand"   # Group log2 fold changes against a grand mean 
#'   type="percentgrand"  # Group log2 fold changes against a grand mean 

#'
#' @param pheno_target List with: 1st object being a character vector
#'   of target column in Pheno 2nd object being a character vector of the target
#'   group(s) in the target Pheno column (1st object).
#'   
#' @param rev Logical whether pairwise comparisons (e.g. log2FC) should be
#'   reversed (default=FALSE).
#'   
#' @param merge_pac Logical whether simplified annotation vector should
#'   automatically be added to the Anno object of the input PAC list object
#'   (default=TRUE). If \code{merge_pac=FALSE} a data.frame is returned. 
#'   
#' @return A PAC object with a pheno_summary folder containing the summarized
#'   data in a data.frame. The data.frame will be named according to the
#'   pheno_target, type and norm input.
#'
#' @examples
#' 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                  package = "seqpac", mustWork = TRUE))
#' 
#' PAC_check(pac) # TRUE
#' 
#' # Easy to generate simple group summaries 
#' pac <- PAC_summary(pac, norm = "cpm", 
#'                    type = "means", pheno_target=list("stage"))       
#' pac <- PAC_summary(pac, norm = "cpm", 
#'                    type = "se", pheno_target=list("stage"))
#' pac <- PAC_summary(pac, norm = "cpm", 
#'                    type = "log2FC", pheno_target=list("stage"))
#' 
#' names(summary(pac))               # Names of individual summaries
#' head(summary(pac)$cpmMeans_stage) # View individual summaries
#' summary(pac)  # View merge summaries
#' df <- as.data.frame(tibble::as_tibble(summary(pac))) # Merge multiple summaries
#' head(df)
#' 
#' 
#' # If a pheno_target is left out, a mean of all samples will be returned:
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                  package = "seqpac", mustWork = TRUE))
#' pac <- PAC_summary(pac, norm = "cpm", type = "mean")  
#' pac <- PAC_summary(pac, norm = "cpm", type = "percentgrand")
#' names(summary(pac))
#' summary(pac)   
#' 
#' @export
#' 
PAC_summary <- function(PAC, norm="counts", type="means", pheno_target=NULL, 
                        rev=FALSE, merge_pac=TRUE){
  
## Prepare and subset ################
  ## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  
  if(!is.null(pheno_target)){ 
    if(length(pheno_target)==1){ 
      pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))
    }
  }
  PAC -> sav
  PAC <- PAC_filter(PAC, subset_only=TRUE, 
                    pheno_target=pheno_target, anno_target=NULL)

  ### Extract data ###
  if(norm=="counts"){
    data <- PAC$Counts
  }else{  
    if(is.null(PAC$norm[[norm]])){
      stop("\nThere is no object called '", norm, "' in the norm list.",
           "\n(Hint: Did you forget to normalize the data using for",
           "\nexample PAC_rpm, or would you rather run the function",
           "\non raw counts using norm='counts'?)")
      }  
    data <- PAC$norm[[norm]]
  }
  
  ### Subset dataset ###
  pheno <- PAC$Pheno
  pheno$All <- "All"
  if(is.null(pheno_target)){
    warning("No grouping factor was specified with pheno_target.",
            "\nCalculations are based on all samples.")
    pheno_target <- list("All","All")
  }else{
    if(length(pheno_target)==1){
      pheno_target[[2]] <- unique(pheno[, pheno_target[[1]]])
    }else{
      if(is.null(pheno_target[[2]])){
        pheno_target[[2]] <- unique(pheno[, pheno_target[[1]]])
      }}}
  indx <- pheno[, pheno_target[[1]]] %in% pheno_target[[2]]
  pheno <- pheno[indx,,drop=FALSE]
  data <- data[,indx,drop=FALSE]

  ### If a group has only one entrance double ###
  ph <- pheno[, pheno_target[[1]]]
  ph <- ph[ph %in% pheno_target[[2]]]
  o_one<- table(ph) ==1
  if(any(o_one)) {
  one_nam  <- names(o_one[o_one==TRUE])
  warning("Only found one sample in at least one group",
          "\nSummary is based on this sample only!")
  dup_dat <- data[,ph %in% one_nam, drop=FALSE]
  samp_one <- colnames(data)[ph %in% one_nam]
  new_nam <-  paste0(samp_one, "(dup)")
  colnames(dup_dat) <- new_nam
  data <- cbind(data, dup_dat)
  dup_ph<- pheno[pheno$Sample_ID %in%  samp_one,, drop=FALSE]
  rownames(dup_ph) <- new_nam
  dup_ph$Sample_ID <- new_nam
  pheno<- rbind(pheno, dup_ph)
  }
  #Fix list with goup data
  stopifnot(identical(rownames(pheno), colnames(data)))
  start_lst <- as.list(pheno_target[[2]])
  names(start_lst) <- pheno_target[[2]]
  sub_data_lst <- lapply(start_lst, 
                         function(x){ 
                           y <- data[, pheno[, pheno_target[[1]]] == x]
                           return(y)
                         })
  ### Create pairwise combinations of pheno_target###
  if (length(sub_data_lst) > 1) { 
    combn_lst <- as.list(data.frame(utils::combn(seq.int(length(sub_data_lst)),m = 2)))
    }
  if (rev == TRUE) { combn_lst <- lapply(combn_lst, function(x) { 
    c(x[2], x[1]) })
  }
  
  ### Apply log2FC to all pairwise combinations ###
  if(type %in% c("log2FC", "Log2FC", "Log2Fc")){
    group_means <- lapply(sub_data_lst, function(x){ 
      mns <- rowMeans(x)
      mns[mns == 0] <- 0.0000001
      return(mns)
      })
    log2FC_lst <- lapply(combn_lst, function(x){
      combn_data <- do.call("cbind", group_means[x])
      log2FC <- as.data.frame(log2(combn_data[,1]/combn_data[,2]))
      colnames(log2FC) <- paste(names(group_means)[x], collapse="_vs_")
      return(log2FC)
      })
    fin <- do.call("cbind", log2FC_lst)
    rownames(fin) <- rownames(data)
    df_nam <- paste0("Log2FC_", pheno_target[[1]])
  }
  if(type %in% c("log2FCgrand", "Log2FCgrand", "log2FCGrand")){
    grnd_mns <- rowMeans(data)
    grnd_mns[grnd_mns == 0] <- 0.0000001
    group_means <- lapply(sub_data_lst, function(x){ mns <- rowMeans(x)
    mns[mns == 0] <- 0.0000001
    return(mns)})
    log2FC_lst <- lapply(group_means, function(x){
      combn_data <- data.frame(group_means=x, grand_means=grnd_mns)
      log2FC <- as.data.frame(log2(combn_data[,1]/combn_data[,2]))
      return(log2FC)
    })
    fin <- do.call("cbind", log2FC_lst)
    rownames(fin) <- rownames(data)
    colnames(fin) <- paste0(names(log2FC_lst), "_vs_grandMeans")
    df_nam <- paste0("Log2FCgrand_", pheno_target[[1]])
  }
  
  ### Apply perc_change to all pairwise combinations ###       
  if(type %in% c("percent", "perc")){
    group_means <- lapply(sub_data_lst, function(x){ 
      as.data.frame(rowMeans(x))
      })
    perc_diff_lst <- lapply(combn_lst, function(x){
      combn_data <- do.call("cbind", group_means[x])
      perc_diff <- as.data.frame((combn_data[,1]/combn_data[,2])*100)
      colnames(perc_diff) <- paste(names(group_means)[x], collapse="_vs_")
      return(perc_diff)
    })
    fin <- do.call("cbind", perc_diff_lst)
    rownames(fin) <- rownames(data)
    df_nam <- paste0("perc_", pheno_target[[1]])
  }
  if(type %in% c("percentgrand", "percgrand", "percentGrand", "percGrand")){
    grnd_mns <- rowMeans(data)
    group_means <- lapply(sub_data_lst, function(x){ 
      as.data.frame(rowMeans(x))
      })
    perc_diff_lst <- lapply(group_means, function(x){
      combn_data <- data.frame(group_means=x, grand_means=grnd_mns)
      perc_diff <- as.data.frame((combn_data[,1]/combn_data[,2])*100)
      return(perc_diff)
    })
    fin <- do.call("cbind", perc_diff_lst)
    rownames(fin) <- rownames(data)
    colnames(fin) <- paste0(names(perc_diff_lst), "_vs_grandMeans")
    df_nam <- paste0("percGrand_", pheno_target[[1]])
  }
  
  
### means, sd and se ###  
  if(type %in% c("means", "Means", "mean", "Mean")){
    group_means <- lapply(sub_data_lst, function(x){ 
      as.data.frame(rowMeans(x))
      })
    fin <- do.call("cbind", group_means)
    colnames(fin) <- names(group_means)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"Means_", pheno_target[[1]])
  }
  if(type %in% c("sd", "SD")){
    group_sd <- lapply(sub_data_lst, function(x){
      as.data.frame(apply(x, 1, stats::sd))
      })
    fin <- do.call("cbind", group_sd)
    colnames(fin) <- names(group_sd)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"SD_", pheno_target[[1]])
  }
  if(type %in% c("se", "SE")){
    group_se <- lapply(sub_data_lst, function(x){
      as.data.frame(apply(x, 1, function(y){stats::sd(y)/sqrt(length(y))}))
      }) 
    fin <- do.call("cbind", group_se)
    colnames(fin) <- names(group_se)
    rownames(fin) <- rownames(data)
    df_nam <- paste0(norm,"SE_", pheno_target[[1]])
  }     
  
  ### Add more summary functions###
##################################################
  
  ### Fix names and return object t###
  if(merge_pac==TRUE){
    PAC <- sav
    PAC$summary$new <- list(NULL)
    PAC$summary[[which(names(PAC$summary)=="new")]] <- as.data.frame(fin)
    
    if(df_nam %in% names(PAC$summary)){
      if(rev==TRUE){
        df_nam<-stringr::str_c(df_nam,"_rev")
       }else{
      numb<-(names(PAC$summary) %in% df_nam)
      numb2<-length(numb[numb==TRUE])
      df_nam<-stringr::str_c(df_nam,"_",numb2)
      }
    }
    
    names(PAC$summary)[which(names(PAC$summary)=="new")] <- df_nam
    stopifnot(PAC_check(PAC))
    if(tp=="S4"){
      return(as.PAC(PAC))
    }else{
      return(PAC)
    }  
  }else{
    fin_lst <- list(fin)
    names(fin_lst) <- df_nam
    return(fin_lst)
  }
}



  
