#' Filter a PAC object on sequence size and coverage  
#'
#' \code{PAC_filter} Filter PAC objects.
#'
#' Given a PAC object the function will extract sequences within a given size
#' interval and percent coverage across independent samples.
#' 
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a Counts table with raw counts or counts per million (cpm).
#'   
#' @param size Integer vector giving the size interval, as c(min,max), that
#'   should be saved (default=c(min,max)).
#'   
#' @param threshold Integer giving the threshold in counts or normalized counts
#'   that needs to be reached for a sequence to be included (default=0).
#'   
#' @param coverage Integer giving the percent of independent samples that need
#'   to reach the threshold for a sequence to be included (default=0).
#'   
#' @param norm Character specifying if filtering should be done using "counts",
#'   "cpm" or another normalized data table in PAC$norm (default="counts").
#' 
#' @param stat (optional) Logical specifying if a coverage graph should be
#'   generated and if users should be prompted prior to proceeding.
#'   (default=FALSE).
#'
#' @param subset_only Logical whether only subsetting using pheno_target and/or
#'   anno_target should be done. If subset=FALSE (default) both subsetting and
#'   other filtering will be done.
#'   
#' @param pheno_target (optional) List with: 
#'          1st object being a character vector of target column in Pheno, 
#'          2nd object being a character vector of the target group(s) in 
#'          the target Pheno column (1st object).
#'          (default=NULL)
#'
#' @param anno_target (optional) List with: 
#'          1st object being a character vector of target column in Anno, 2nd
#'          object being a character vector of the target type/biotypes(s) in
#'          the target Anno column (1st object).
#'          (default=NULL) 
#'          
#' @param allbut (optional) TRUE or FALSE to give option to include
#'          all possible variables but the ones defined in pheno_target
#'          and anno_target. Useful when you wish to for example filter out
#'          one particular variable.
#'
#' @return A list of objects: 
#'               PAC object with filtered data.   
#'               (optional) A coverage plot 
#'               
#' @examples
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                  package = "seqpac", mustWork = TRUE))
#' 
#'###--------------------------------------------------------------------- 
#'## Extracts all sequences between 10-80 nt in length with at 
#'## least 5 counts in 20% of all samples.
#'pac_lowfilt <- PAC_filter(pac, size=c(10,80), threshold=5,
#'                          coverage=20, norm = "counts",
#'                          pheno_target=NULL, anno_target=NULL)
#'
#'###--------------------------------------------------------------------- 
#'## Extracts sequences with 22 nt size and the samples in Batch1 and Batch2.
#'pac_subset <- PAC_filter(pac, subset_only = TRUE,
#'                         pheno_target=list("batch", c("Batch1", "Batch2")), 
#'                         anno_target=list("Size", "22"))
#'
#'###--------------------------------------------------------------------- 
#'## Extracts all sequences with >=5 counts in 100% of samples a within stage
#'filtsep <- PAC_filtsep(pac, norm="counts", threshold=5, 
#'                       coverage=100, pheno_target= list("stage"))
#'
#'pac_filt <- PAC_filter(pac, subset_only = TRUE,
#'                      anno_target= unique(do.call("c", as.list(filtsep))))
#'                      
#'###--------------------------------------------------------------------- 
#'## Extract all biotypes but tRNA in all samples but in sample fastq_1
#'
#'pac_filt <- PAC_filter(pac, anno_target=list("Biotypes_mis0","tRNA"),
#'                      pheno_target=list("sample","fastq1"), allbut=TRUE)
#' 
#' 
#' 
#' @export

PAC_filter <- function(PAC, size=NULL, threshold=0, coverage=0, 
                       norm="counts", subset_only=FALSE, stat=FALSE, 
                       pheno_target=NULL, anno_target=NULL, allbut=FALSE){
  

  ## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  
  options(scipen=999)
  
  
  strt <- nrow(PAC$Counts)
  nsamp <- ncol(PAC$Counts)
  x_graph <- n_features <- NULL
  
  if(allbut==TRUE){
    pheno_t<-(unique(PAC$Pheno[,pheno_target[[1]]]))
    pheno_target[[2]] <- pheno_t[pheno_t != pheno_target[[2]]]
    anno_t <-(unique(PAC$Anno[,anno_target[[1]]]))
    anno_target[[2]]<-anno_t[anno_t != anno_target[[2]]]
  }
  
  ### Subset samples by Pheno 
  if(!is.null(pheno_target)){
    if(length(pheno_target)==1){
      pheno_target[[2]] <- unique(PAC$Pheno[, pheno_target[[1]]])
    }
    sub_pheno <- as.character(
      PAC$Pheno[, pheno_target[[1]]]) %in% pheno_target[[2]]
    if(any(names(PAC)=="norm")){
      PAC$norm <- lapply(as.list(PAC$norm), function(x){
        x[,sub_pheno, drop=FALSE]
      })
    }
    if(any(names(PAC)=="summary")){
      if(!any(sub_pheno)){
        warning(
          "\nTable(s) were found in PAC@summary that may have been generated",
          "\nwith samples that now are removed. Summary table names has ",
          "\ntherefore been marked.")
        names(PAC$summary) <- paste0(names(PAC$summary), 
                                     "_WARNING_pheno_filter")
      }}
    PAC$Counts  <- PAC$Counts[,sub_pheno, drop=FALSE]
    PAC$Pheno  <- PAC$Pheno[sub_pheno,, drop=FALSE]
    ord_mch <- order(match(
      as.character(PAC$Pheno[, pheno_target[[1]]]), pheno_target[[2]]))
    PAC$Counts <- PAC$Counts[,ord_mch,drop=FALSE]
    PAC$Pheno <- PAC$Pheno[ord_mch,,drop=FALSE]
    if(any(names(PAC)=="norm")){
      PAC$norm <- lapply(PAC$norm, function(x){x[, ord_mch, drop=FALSE]})
    }
    tab_pheno <- as.data.frame(table(sub_pheno))
    passed <- tab_pheno[tab_pheno[,1]==TRUE, 2]
    if(passed==0){
      stop("Pheno filter resulted in 0 samples.")
    }else{
      cat(paste0("\n-- Pheno target was specified, will retain: ", 
                 passed, " of ", length(sub_pheno), " samples."))
    }
  }  
  
  ### Subset data by Anno
  if(!is.null(anno_target)){
    if(methods::is(anno_target,"list")){
      if(length(anno_target)==1){
      anno_target[[2]] <- unique(PAC$Anno[, anno_target[[1]]])
     }
    }else{
      PAC$Anno$all_names <- rownames(PAC$Anno)
      anno_target <- list("all_names", anno_target)
    }
    sub_anno <- as.character(
      PAC$Anno[,anno_target[[1]]]) %in% as.character(anno_target[[2]])
    if(any(names(PAC)=="norm")){
      PAC$norm <- lapply(as.list(PAC$norm), function(x){
        x[sub_anno,, drop=FALSE]})}
    if(any(names(PAC)=="summary")){
      PAC$summary <- lapply(as.list(PAC$summary), function(x){
        x[sub_anno,, drop=FALSE]})}
    PAC$Counts  <- PAC$Counts[sub_anno,,drop=FALSE]
    PAC$Anno  <- PAC$Anno[sub_anno,,drop=FALSE]
    
    ord_mch <- order(match(
      as.character(PAC$Anno[,anno_target[[1]]]), anno_target[[2]]))
    PAC$Counts <- PAC$Counts[ord_mch,,drop=FALSE]
    PAC$Anno <- PAC$Anno[ord_mch,,drop=FALSE]
    if(any(names(PAC)=="norm")){
      PAC$norm <- lapply(PAC$norm, function(x){
        x[ord_mch,, drop=FALSE]
      })
    }
    if(any(names(PAC)=="summary")){
      PAC$summary <- lapply(PAC$summary, function(x){
        x[ord_mch,, drop=FALSE]
      })
    }
    tab_anno <- as.data.frame(table(sub_anno))
    passed <- tab_anno[tab_anno[,1]==TRUE, 2]
    if(passed==0){
      stop("Anno filter resulted in 0 sequences.") 
    }else{
      cat(paste0("\n-- Anno target was specified, will retain: ", 
                 passed, " of ", length(sub_anno), " seqs."))
    }
  }  
  
  ### Subset data by Size
  if(!is.null(size)){
    if(!length(size)==2)
      stop("\nYou must specify both min and max size.",
           "\nOnly one number was obtained.")
    sub_size <- PAC$Anno$Size >= size[1] & PAC$Anno$Size <= size[2]
    if(any(names(PAC)=="norm")){
      PAC$norm <- lapply(as.list(PAC$norm), function(x){
        x[sub_size, , drop=FALSE]
      })
    }
    if(any(names(PAC)=="summary")){
      PAC$summary <- lapply(as.list(PAC$summary), function(x){
        x[sub_size, , drop=FALSE]
      })
    }
    PAC$Counts  <- PAC$Counts[sub_size, , drop=FALSE] 
    PAC$Anno  <- PAC$Anno[sub_size, , drop=FALSE]
    tab_anno <- as.data.frame(table(sub_size))
    passed <- tab_anno[tab_anno[,1]==TRUE, 2]
    if(passed==0){
      stop("Size filter resulted in 0 sequences.") 
    }else{
      cat(paste0("\n-- Size filter will retain: ", 
                 passed, " of ", length(sub_size), " seqs."))
    }
  }
  if(!subset_only==TRUE){   
    ### Extract essential information
    if(!norm %in% c(names(PAC$norm),"counts", "Counts")){
      stop("\nThe data specified in 'norm' was not avaiable in PAC.",
           "\nMake sure that name in 'norm' specifies a table name in",
           "\nPAC@Counts or PAC@norm$...")
    }
    if(norm %in% c("counts","Counts")){ 
      df <- PAC$Counts
      cat("\n-- Count filter was specified.")
    }else{ 
      df <- PAC$norm[[norm]]
      cat(paste0("\n-- Filter on normalized (", 
                 norm, ") values was specified."))
    }   
    
    ### Check col and row names
    if(!identical(rownames(PAC$Anno), rownames(df))){
      stop("\nSequence names not matching in input files! (Anno vs Counts)")
    }
    if(!identical(rownames(PAC$Pheno), colnames(df))){
      stop("\nSample names not matching in input files!(Pheno vs Counts)")
    }
    
    ### Create filter
    idx_filt <- data.frame(
      rowSums(df >= threshold)) >= round(ncol(df)*(coverage*0.01))
    idx_tab <- as.data.frame(table(idx_filt))
    
    ### Generate stat
    if(stat==TRUE){
      filt_plot <- data.frame(matrix(NA, ncol=2, nrow=9))
      colnames(filt_plot) <- c("filter","n_features")
      filt_plot$x_graph <- as.numeric(c(1,5,10,20,30,40,50,75,100))
      filt_plot$filter <- c(paste0(">1_in_", coverage, "%"), 
                            paste0(">5_in_", coverage, "%"),
                            paste0(">10_in_", coverage, "%"),
                            paste0(">20_in_", coverage, "%"), 
                            paste0(">30_in_", coverage, "%"), 
                            paste0(">40_in_", coverage, "%"), 
                            paste0(">50_in_", coverage, "%"),
                            paste0(">75_in_", coverage, "%"), 
                            paste0(">100_in_", coverage, "%"))
      for (i in seq.int(length(filt_plot$x_graph))){ 
        tab <- as.data.frame(table(data.frame(
          rowSums(df >= filt_plot$x_graph[i])) >= round(
            ncol(df)*(coverage*0.01))))
        filt_plot[i,2] <- tab$Freq[tab$Var1=="TRUE"]
      }
      ### Plot graph
      try_err<- try( 
        p <- ggplot2::ggplot(
          filt_plot, ggplot2::aes(x=x_graph, y=n_features, fill=x_graph))+
          ggplot2::geom_line()+
          ggplot2::geom_point(cex=2, fill="blue")+
          ggplot2::geom_hline(yintercept=0)+
          ggplot2::scale_x_discrete(limit=filt_plot$x_graph, 
                                    labels= as.character(filt_plot$filter))+
          ggplot2::ggtitle("User filter:")+
          ggplot2::xlab(NULL)+
          ggplot2::theme_classic()+
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                             hjust = 0))+
          ggplot2::geom_vline(xintercept=threshold, col="red")+
          ggplot2::geom_label(
                   x=50, y=max(filt_plot[,2])*0.95, 
                   label=paste0("n input sequences: ", strt, 
                         "\nn analyzed sequences: ", nrow(df), 
                         "\nn sequences after filter: ",
                         idx_tab[idx_tab[,1]==TRUE, 2]), show.legend = FALSE)+
          ggplot2::theme(plot.margin=ggplot2::margin(t = 1, r = 1,
                                                     b = 1, l = 1,
                                                     unit="cm"),
                         plot.title = ggplot2::element_text(color="red", 
                                                            size=10)))
      
      if(is(try_err, "try-error")){
        warning("Was unable to create filtering graph. Probable reason:", 
                "\nCorrupt PAC or very few reads from start.")
      }
      
      Sys.sleep(0.01)
      print(p)
      ### Promt for user input
      cat("\n!!           !!\nUser input needed:\n") 
      answer <- readline(prompt="Continute with this filter? [Y/n]")
      if(answer=="Y" | answer=="y"){}else{stop("Script was terminated by user.")
      }
    }
    
    ### Apply filter
    passed <- idx_tab[idx_tab[,1]==TRUE, 2]
    if(passed==0){
      stop("The chosen filter resulted in 0 sequences.") 
    }else{
      cat(paste0("\n-- The chosen filters will retain: ", 
                 idx_tab[idx_tab[,1]==TRUE, 2], " of ", strt, " seqs."))   
      if(any(names(PAC)=="norm")){
        PAC$norm <- lapply(as.list(PAC$norm), function(x){
          x[idx_filt,, drop=FALSE]
          })
      }
      if(any(names(PAC)=="summary")){
        PAC$summary <- lapply(as.list(PAC$summary), function(x){
          x[idx_filt,, drop=FALSE]
          })
      }
      PAC$Counts  <- PAC$Counts[idx_filt,, drop=FALSE] 
      PAC$Anno  <- PAC$Anno[idx_filt,, drop=FALSE]
    }
  }
  ## Double check and return
  
  if(PAC_check(PAC)==TRUE){
    if(tp=="S4"){
       return(as.PAC(PAC))
    }else{
       return(PAC)
    }
  }
}      

