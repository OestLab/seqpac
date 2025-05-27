#' Generates size distribution plots from a PAC object
#'
#' \code{PAC_sizedist} plotting size distribution with bar charts, allowing for
#' visualization of sequence classes and summaries.
#'
#' Given a PAC object the function will attempt to order sequences by their 
#' size (number of nucleotides) and visualize the contribution of specific 
#' classes of sequences (e.g. sRNA classes) at each size point.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a count table with raw counts.
#'
#' @param range Integer vector giving the range  in sequence lengths
#'   (default=c(min, max)).
#'   
#' @param colors Character vector with RGB color codes to be parsed to ggplot2. 
#'
#' @param norm Character indicating what type of data to be used. If
#'  type="counts" the plots will be based on the raw Counts. If type="cpm" the
#'  analysis will be done on cpm values returned from \code{PAC_norm} function
#'  and stored in the norm folder of the PAC-list object. The name of any other
#'  table in the norm(PAC) folder can also be used.  
#' 
#' @param ymax Integer that sets the maximum y for all plots (all plots gets the
#'   same y-axes). If ymax=NULL, then ggplot2 will automatically set ymax for
#'   each plot individually (different y-axes).
#' 
#' @param anno_target List with: 
#'                          1st object being character vector of target
#'                          column(s) in Anno, 2nd object being a character
#'                          vector of the target biotypes(s) in the target
#'                          column (1st object). (default=NULL)
#' 
#' @param pheno_target List with: 
#'                          1st object being character vector of target
#'                          column(s) in Pheno, 2nd object being a character
#'                          vector of the target group(s) in the target column
#'                          (1st object). (default=NULL)
#' 
#' @param summary_target List with: 
#'                          1st object being character target object in
#'                          summary(PAC), 2nd object being a character vector of
#'                          the target columns(s) in the target object (1st
#'                          object). (default=NULL)
#' 
#'
#' @return A list of objects: 
#'               1st object (Histograms::Samples): Individual histograms showing
#'               the nucleotide ratios per sample over the specified range.
#'               2nd object (Data::Samples): Data used to generate the plots.
#'               3rd object (Stacked_bars::Groups): Stacked bars showing the
#'               mean ratios of each nucleotide per group over the specified
#'               range.
#'               4th object (Error_bars::Groups): Error bar plots with mean
#'               ratio of each nucleotide per group over the specified range.
#'
#' @examples
#' 
#' ##########################################
#' ### Stacked bars in seqpac 
#' ##----------------------------------------
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' 
#' PAC_filt <- PAC_norm(pac, norm="cpm")
#' PAC_filt <- PAC_summary(PAC=PAC_filt, norm = "cpm", 
#'                         type = "means", pheno_target=list("stage"))
#' 
#' 
#' ord <- c("no_anno", "other", "miRNA", "tRNA", "rRNA", "snoRNA",  "lncRNA")
#' 
#' sizedist_plots <- PAC_sizedist(PAC_filt, 
#'                                anno_target=list("Biotypes_mis0", ord), 
#'                                summary_target=list("cpmMeans_stage"))
#' cowplot::plot_grid(plotlist=sizedist_plots[[1]], nrow = 3, ncol = 1)
#' 
#' sizedist_plots <- PAC_sizedist(PAC_filt, norm="counts", 
#'                                anno_target=list("Biotypes_mis0", ord), 
#'                                pheno_target=list("batch", "Batch1"))
#' cowplot::plot_grid(plotlist=sizedist_plots[[1]], nrow = 3, ncol = 1)
#' 
#' 
#' @export

PAC_sizedist <- function(PAC, norm="counts", range=NULL, anno_target, 
                         pheno_target=NULL, summary_target=NULL, colors=NULL,
                         ymax=NULL){
  
  
  size <- biotype <- NULL
  # Prepare filtered PAC
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
  
  if(!is.null(anno_target)){ 
    if(length(anno_target)==1){
      anno_target[[2]] <- as.character(unique(PAC$Anno[,anno_target[[1]]]))
     }
    }
  if(is.null(range)){
    range <- c(min(PAC$Anno$Size), max(PAC$Anno$Size))
    }
  
  PAC <- seqpac::PAC_filter(PAC, size=range, pheno_target=pheno_target, 
                            anno_target=anno_target, subset_only=TRUE)
  
  stopifnot(PAC_check(PAC))
  ph <- PAC$Pheno
  anno <- PAC$Anno
  
  # Extract data
  if(is.null(summary_target)){
    if(is.null(norm)){
      data <- PAC$Counts; labl <- "Counts"
    }else{
      if(norm=="counts"){
        data <- PAC$Counts; labl <- "Counts"
      }else{
        data <- PAC$norm[[norm]]; labl <- norm
      }
    }
  }else{
    data <- PAC$summary[[summary_target[[1]]]]; labl <- summary_target
    }
  

  #### Summarize over size and biotype
  bio_fact <- factor(anno[, anno_target[[1]]], levels=anno_target[[2]])
  seq_range <- seq(range[1], range[2])
  size_fact <- factor(anno$Size, levels=seq_range)
  size_lst <- lapply(as.list(data), function(x){
    bio_agg <- stats::aggregate(x, list(bio_fact, size_fact), sum)
    colnames(bio_agg) <- c("biotype", "size", "data")
    bio_agg_lst <- lapply(split(bio_agg, bio_agg$biotype), function(y){
      if(any(!seq_range %in% y$size )){
        y <- rbind(y, 
                   data.frame(biotype=as.character(unique(y$biotype)), 
                              size=seq_range[
                                !seq_range %in% as.character(y$size)], 
                              data=0))
      }
      fin <- y[order(y$size),,drop=FALSE]
      stopifnot(identical(as.character(fin$size), as.character(seq_range)))
      return(fin)
    })
    fin_df <- do.call("rbind", bio_agg_lst)
    rownames(fin_df) <- NULL
    return(fin_df)
  })
  
  #### Generate individial histograms 
  options(scipen=999)
  #### Set up colors colors ###
  if(is.null(colors)){
    colfunc <- grDevices::colorRampPalette(c("#094A6B", "#EBEBA6", "#9D0014"))
    rgb <- colfunc(sum(!anno_target[[2]] %in% c("other", "no_anno")))
    rgb_vec <- NULL
    cnt <- 0
    for(i in seq.int(length(anno_target[[2]]))){ 
      if(anno_target[[2]][i] == "other"){
        rgb_vec[i] <- "#808080";  cnt <- cnt+1}
      if(anno_target[[2]][i] == "no_anno"){
        rgb_vec[i] <- "#C0C0C0";  cnt <- cnt+1}
      if(!anno_target[[2]][i] %in% c("other","no_anno")){
        rgb_vec[i] <- rgb[i-cnt]}
    }
  }else{
    rgb_vec <- colors
    }
  
  #### Plot individual plots ###
  histo_lst <- list(NA)
  if(!is.null(summary_target)){
    samp <- colnames(data)
  }else{
    if(is.null(pheno_target)){
      samp <- rownames(ph)
    }else{
      samp <- paste0(ph[,pheno_target[[1]]],"-", rownames(ph)) 
    }
  }
  for(i in seq.int(length(size_lst))){
    histo_lst[[i]] <- ggplot2::ggplot(size_lst[[i]], 
                                      ggplot2::aes(x=size, 
                                                   y=data, 
                                                   fill=biotype))+
      ggplot2::geom_bar(width = 0.9, size =0.2, colour="black", stat="identity")+
      ggplot2::geom_hline(yintercept=0, col="azure4")+
      
      ggplot2::xlab("Size (nt)")+
      ggplot2::ylab(paste0(labl))+
      ggplot2::labs(subtitle = samp[i])+
      ggplot2::scale_fill_manual(values=rgb_vec)+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
     if(!is.null(ymax)){
         histo_lst[[i]] <- histo_lst[[i]]+ 
           ggplot2::coord_cartesian(ylim=c(0, ymax))
     }
    
    names(histo_lst)[i] <- names(size_lst)[i]
  }
  return(list(Histograms=histo_lst, Data=size_lst))
}

