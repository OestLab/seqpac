PAC_stackbar <- function(PAC, anno_target=NULL, pheno_target=NULL, summary_target=NULL, 
                         color=NULL, width=1.0, no_anno=TRUE, total=TRUE, 
                         norm="counts", style="percent"){
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  
  stopifnot(PAC_check(PAC))
  Sample <- Value <- Category <- tot_counts <- NULL
  
  if(!is.null(pheno_target)){ 
  if(length(pheno_target)==1){ 
  if(is(PAC$Pheno[,pheno_target[[1]]], "factor")){pheno_target[[2]] <- levels(PAC$Pheno[,pheno_target[[1]]])}
  else{pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))}
    }}
  else{
    PAC$Pheno$eXtra_Col <- rownames(PAC$Pheno)
    pheno_target <- list(NA)
    pheno_target[[1]] <- "eXtra_Col"
    pheno_target[[2]] <- PAC$Pheno$eXtra_Col}
  
  if(!is.null(anno_target)){ 
  if(length(anno_target)==1){
  if(is(PAC$Anno[,anno_target[[1]]], "factor")){anno_target[[2]] <- levels(PAC$Anno[,anno_target[[1]]])}
    else{anno_target[[2]] <- as.character(unique(PAC$Anno[,anno_target[[1]]]))}}
  }
  
  PAC_sub <- seqpac::PAC_filter(
    PAC,
    subset_only=TRUE, 
    pheno_target=pheno_target,
    anno_target=anno_target
  )
  
  cat("\n\n")
  
  if(isS4(PAC_sub)){
    PAC_sub <- as(PAC_sub, "list")
  }
  
  anno <- PAC_sub$Anno
  pheno <- PAC_sub$Pheno
  
  if(norm == "counts"){data <- PAC_sub$Counts}
  else{data <- PAC_sub$norm[norm][[1]]}
  
  if(!is.null(summary_target)){data <- PAC_sub$summary[[summary_target[[1]]]]
  if(length(summary_target) > 1){summary_groups <- summary_target[[2]]
  if(!all(summary_groups %in% colnames(data))){missing_groups <- summary_groups[ !summary_groups %in% colnames(data)]
  stop("The following summary_target groups are not present in the summary table: ",paste(missing_groups, collapse=", "),
       "\nAvailable groups are: ",paste(colnames(data), collapse=", "))
      }
      
      data <- data[, summary_groups, drop=FALSE]
    }
  }
  
  if(no_anno==FALSE){
    data <- data[!as.character(anno[, anno_target[[1]]]) == "no_anno",]
    anno <- anno[!as.character(anno[, anno_target[[1]]]) == "no_anno",]
  }
  
  data_shrt <- stats::aggregate(
    data,
    list(anno[, anno_target[[1]]]),
    "sum"
  )
  
  tot_cnts <- colSums(data)
  data_shrt_total <- data_shrt
  data_long_tot <- reshape2::melt(
    data_shrt_total,
    id.vars="Group.1"
  )
  
  colnames(data_long_tot) <- c("Category", "Sample", "Value")
  
  if(style=="percent"){
    data_shrt_perc <- data_shrt
    data_shrt_perc[,-1] <- "NA"
    
    for(i in seq.int(length(tot_cnts))){ 
      data_shrt_perc[,1+i] <- data_shrt[,1+i]/tot_cnts[i]
    }
    
    data_long_perc <- reshape2::melt(
      data_shrt_perc,
      id.vars="Group.1"
    )
    
    colnames(data_long_perc) <- c("Category", "Sample", "Value")
    data_long_perc$Value <- data_long_perc$Value * 100
    percentage_table <- reshape2::dcast(
      data_long_perc,
      Category ~ Sample,
      value.var="Value"
    )
    
    percentage_table[-1] <- lapply(
      percentage_table[-1],
      function(x) paste0(round(x, 2), "%")
    )
    
    print(percentage_table)
    
    data_long_perc <- data_long_perc
    }
  else{
    data_long_perc <- data_long_tot
  }
  
  bio <- anno_target[[2]] 
  extra <- which(bio %in% c("no_anno", "other"))
  
  if(length(extra)>0){
    bio <- c(sort(bio[extra]), bio[-extra])
  }
  
  data_long_perc$Category <- factor(
    as.character(data_long_perc$Category), 
    levels=bio
  )
  
  if(is.null(pheno_target) && is.null(summary_target)){
    data_long_perc$Sample <- factor(
      as.character(data_long_perc$Sample), 
      levels=as.character(unique(data_long_perc$Sample))
    )
  }
  
  tot_cnts <- tot_cnts[
    match(names(tot_cnts), unique(data_long_perc$Sample))
  ]
  
  data_long_perc$tot_counts <- ""
  
  if(total==TRUE){
    trg_1st <- levels(data_long_perc$Category)[
      length(levels(data_long_perc$Category))
    ]
    
    data_long_perc$tot_counts[
      data_long_perc$Category == trg_1st
    ] <- tot_cnts
  }
  
  if(is.null(color)){
    n_extra <- length(extra)
    colfunc <- grDevices::colorRampPalette(
      c("#094A6B", "#EBEBA6", "#9D0014")
    )
    
    if(n_extra==1){
      color <- c(colfunc(length(bio)-1), "#6E6E6E")
    }
    
    if(n_extra==2){
      color <- c(
        colfunc(length(bio)-2),
        "#6E6E6E",
        "#BCBCBD"
      )
    }
    
  if(n_extra==0){
      color <- colfunc(length(bio))}}
  else{
    color <- rev(color)
  }
  
  p1 <- ggplot2::ggplot(
    data_long_perc,
    ggplot2::aes(x=Sample, y=Value, fill=Category)
  ) +
    ggplot2::geom_bar(
      stat="identity",
      col="black",
      width=width,
      linewidth=0.3
    ) + 
    ggplot2::geom_text(
      ggplot2::aes(label=tot_counts),
      nudge_y=-3,
      nudge_x=0,
      angle=0,
      color="black",
      size=4
    ) +
    ggplot2::geom_hline(
      yintercept=0,
      col="black"
    ) +
    {
      if(style=="percent")
        ggplot2::coord_cartesian(ylim=c(-2, 100))
    } + 
    {
      if(style=="percent")
        ggplot2::ylab("Percent of total reads")
    } +
    {
      if(style=="total")
        ggplot2::ylab("Total reads")
    } +
    {
      if(style=="percent")
        ggplot2::geom_hline(
          yintercept=100,
          col="black"
        )
    } +
    ggplot2::scale_fill_manual(
      values=rev(color)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.ticks.length.y=ggplot2::unit(.25, "cm"),
      plot.caption=ggplot2::element_text(size=12, face="bold"),
      axis.title.y=ggplot2::element_text(size=16, face="bold"),
      axis.line=ggplot2::element_blank(),      
      axis.title.x=ggplot2::element_blank(), 
      axis.text=ggplot2::element_text(size=12),
      axis.text.x=ggplot2::element_text(angle=45, hjust=1),
      panel.background=ggplot2::element_blank()
    )
  
  return(p1)
}