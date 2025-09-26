#' Generates a nucleotide bias analysis from a PAC object
#'
#' \code{PAC_nbias} analyses nucleotide bias.
#'
#' Given a PAC object the function will attempt to extract the ratios of
#' specific nucleotides at a given position in sequences in the Anno data.frame
#' in relation to the sequence counts in Counts.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC-list object containing an Anno data.frame with sequences as
#'   row names and a count table with raw counts.
#' @param position Integer indicating the nucleotide position from 3' to 5'
#'   position (default=1).
#' @param nucleotide_range Integer vector indicating the sequence size range
#'   (default=c(min, max)).
#'
#' @param norm Character indicating what type of data to be used. If
#'  type="counts" the plots will be based on the raw Counts. If type="cpm" the
#'  analysis will be done on cpm values returned from \code{PAC_norm} function
#'  and stored in the norm folder of the PAC-list object. The name of any other
#'  table in the norm(PAC) folder can also be used.    
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
#'                          1st object being character vector of target object
#'                          in summary(PAC), 2nd object being a character vector
#'                          of the target column(s) in the target summary object
#'                          (1st object). (default=NULL)
#'
#' @param colors Character vector with RGB color codes to be parsed to ggplot2.
#'
#' @param ymax Integer that sets the maximum y for all plots (all plots gets the
#'   same y-axes). If ymax=NULL, then ggplot2 will automatically set ymax for
#'   each plot individually (different y-axes).
#'   
#' @param data_only logical. If data_only=TRUE a data.frame a simple Anno object
#'   is returned with a Size and a Nucleotide bias column. As default,
#'   data_only=FALSE then graphs are returned in addition to data.
#'   
#' @return A list of objects: 
#'               1st object (Histograms::Samples): Individual histograms showing
#'               the nucleotide ratios per sample over the specified range. 2nd
#'               object (Data::Samples): Data used to generate the plots.
#'               
#' @examples
#' 
#' 
#' # Load a PAC-object 
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                   package = "seqpac", mustWork = TRUE))
#' output_nbias <- PAC_nbias(pac)
#' cowplot::plot_grid(plotlist=output_nbias$Histograms)
#' 
#' # Only miRNA (Oops, heavy T-bias on 1st nt; are they piRNA?)  
#' table(anno(pac)$Biotypes_mis0)
#' output_nbias <- PAC_nbias(pac, anno_target = list("Biotypes_mis0", "miRNA") )
#' cowplot::plot_grid(plotlist=output_nbias$Histograms)
#' 
#' # Switch to 10:th nt bias 
#' output_nbias <- PAC_nbias(pac, position=10, 
#'                           anno_target = list("Biotypes_mis0", "miRNA"))
#' cowplot::plot_grid(plotlist=output_nbias$Histograms)
#' 
#' # Summarized over group cpm means
#' pac_test <- PAC_summary(pac, norm = "cpm", type = "means", 
#'                         pheno_target=list("stage"), merge_pac=TRUE)
#' output_nbias <- PAC_nbias(pac_test, summary_target = list("cpmMeans_stage") )
#' cowplot::plot_grid(plotlist=output_nbias$Histograms)
#' 
#' 
#' @export


PAC_nbias <- function(PAC, position=1, norm=NULL, nucleotide_range=NULL, anno_target=NULL, 
                      pheno_target=NULL, summary_target=NULL, colors=NULL,
                      ymax=NULL, data_only=FALSE){
  
  counts <- nucleotide <- NULL
  
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
  if(is.null(nucleotide_range)){
    nucleotide_range <- c(min(PAC$Anno$Size), max(PAC$Anno$Size))
    }
  
  PAC <- seqpac::PAC_filter(PAC, nucleotide_range=nucleotide_range, pheno_target=pheno_target, 
                            anno_target=anno_target, subset_only=TRUE)
  stopifnot(PAC_check(PAC))
  
  
  # Extract data
  if(is.null(summary_target)){
    if(is.null(norm)){
      dat <- PAC$Counts; labl <- "Counts"
    }
    else{
      if(norm=="counts"){
        dat <- PAC$Counts; labl <- "Counts"
      }
      else{
        dat <- PAC$norm[[norm]]; labl <- norm
      }
    }
  }else{
    dat <- PAC$summary[[summary_target[[1]]]]; labl <- summary_target
    }
  
  # Extract nt anno
  anno <-PAC$Anno
  anno$nuc_bias <- substr(rownames(anno), start=position, stop=position)
  if(data_only==TRUE){
    colmn1  <- which(colnames(anno)=="nuc_bias")
    colmn2  <- which(colnames(anno)=="Size")
    colnames(anno)[colmn1] <- paste0(position, "_nuc_bias") 
    return(anno[,c(colmn2, colmn1),drop=FALSE])
  }else{
   cat(paste0("\nCounting nucleotides"))
   combin <- c(paste0(nucleotide_range[1]:nucleotide_range[2], "_A"),
               paste0(nucleotide_range[1]:nucleotide_range[2], "_T"),
               paste0(nucleotide_range[1]:nucleotide_range[2], "_C"),
               paste0(nucleotide_range[1]:nucleotide_range[2], "_G"),
               paste0(nucleotide_range[1]:nucleotide_range[2], "_N"))
   nuc_lst <- lapply(as.list(dat), function(x){
     nuc_agg <- stats::aggregate(x, list(factor(paste(anno$Size, 
                                               anno$nuc_bias, sep="_"))), sum)
     colnames(nuc_agg) <- c("position_nuc", "counts")
     zeros <- combin[!combin %in% as.character(nuc_agg$position_nuc)]
     if(length(zeros) >0){
       nuc_agg <- rbind(nuc_agg, data.frame(position_nuc=zeros, counts=0))}
     nuc_agg_ord <- nuc_agg[match(combin, nuc_agg$position_nuc),,drop=FALSE]
     rownames(nuc_agg_ord) <- NULL
     stopifnot(identical(as.character(nuc_agg_ord$position_nuc), combin))
     splt<- do.call("rbind", 
                   strsplit(as.character(nuc_agg_ord$position_nuc), split="_"))
    return(data.frame(length=splt[,1], nucleotide=splt[,2], 
                      counts=nuc_agg_ord[,2]))
  })
  
  #### Set options and load requirements
  options(scipen=999)
  #### Set up colors colors ###
  if(is.null(colors)){
    colfunc_sports <- grDevices::colorRampPalette(c("#094A6B", 
                                                    "#FFFFCC", 
                                                    "#9D0014"))
    colors <- colfunc_sports(5)
    colors <- c("#A0A0A0",colors[c(1,2,3,5)])
  }
  
  #### Plot ###                 
  histo_lst <- list(NA)
  if(is.null(summary_target)){samp <- rownames(PAC$Pheno)}
  else{samp <- names(PAC$summary[[summary_target[[1]]]])}
  for(i in seq.int(length(nuc_lst))){
    nuc_lst[[i]]$nucleotide <- factor(nuc_lst[[i]]$nucleotide, 
                                      levels=c("N","C","G","A","T"))
    uni_chr_len <- as.integer(unique(as.character(nuc_lst[[i]]$length)))
    nuc_lst[[i]]$length <- factor(nuc_lst[[i]]$length, 
                                  levels=uni_chr_len[order(uni_chr_len)] )
    histo_lst[[i]] <- ggplot2::ggplot(nuc_lst[[i]], 
                                      ggplot2::aes(x=length, y=counts, 
                                                   fill=nucleotide))+
      ggplot2::geom_bar(width = 0.9, cex=0.2, colour="black", stat="identity")+
      ggplot2::geom_hline(yintercept=0, col="azure4")+
      ggplot2::xlab("Size (nt)")+
      ggplot2::ylab(paste(labl))+
      ggplot2::labs(subtitle = samp[i])+
      ggplot2::scale_fill_manual(values=colors)+
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
    names(histo_lst)[i] <- names(nuc_lst)[i]
    if(!is.null(ymax)){
         histo_lst[[i]] <- histo_lst[[i]]+ 
           ggplot2::coord_cartesian(ylim=c(0, ymax))
     }
  }
  return(list(Histograms=histo_lst, Data=nuc_lst))
  }
}
