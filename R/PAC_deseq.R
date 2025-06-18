#' A wrapper to DESeq2 that apply differential expression analysis on a PAC
#' object
#'
#' \code{PAC_deseq} DESeq2 analysis on PAC_object.
#'
#' Given a PAC object this function will apply a differential expression
#' analysis using DESeq2.
#'
#' @family PAC analysis
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param PAC PAC object containing a Pheno data.frame with samples as row
#'   names, an Anno data.frame with sequences as row names and a Counts table
#'   with raw counts. The Counts table must have the sample names as column
#'   names and sequences as row names.
#'
#' @param model Character vector describing the statistical model based on
#'   column names in Pheno.
#'
#' @param deseq_norm Logical whether to return deseq normalized values or raw
#'   counts (default=FALSE).
#'
#' @param pheno_target (optional) List with: 1st object being a character
#'   indicating the main target column in Pheno. 2nd object being a character
#'   vector of the target group(s) in the target Pheno column (1st object).
#'
#'   Important: In \code{PAC_deseq}, pheno_target controls the main comparative
#'   factor category from which a summarized table and plots will be generated.
#'   If, for instance, a target column named "groups" in pheno(PAC) contains
#'   "control" and "treatment" categories, setting pheno_target=list("groups",
#'   c("treatment", "controls") ensures that "treatment" is presented 1st in the
#'   factor levels, making for example log2FC appear as "treatment vs control".
#'   As default, pheno_target=NULL will result in the factor levels being
#'   automatically sorted in ascending order, which in the example above would
#'   result in control vs treatment. If no pheno_target is given, the
#'   first feature in the model will also be the main comparison presented in
#'   the graphs and summary tables.
#'
#' @param test Character parsed directly to \code{\link[DESeq2]{DESeq}} that
#'   controls what type of statistical test that should be used. Alternatives
#'   are either "Wald" (Wald significance test) or "LTR" (likelihood ratio
#'   test/chi-squared test). See \code{\link[DESeq2]{DESeq}} for more details.
#'   (Default="Wald")
#'
#' @param fitType Character parsed directly to \code{\link[DESeq2]{DESeq}} that
#'   controls what type of dispersion fit that should be used. Alternatives are
#'   either "parametric" (dispersion-mean relation), "local" (local regression
#'   of log dispersions), "mean" (mean of gene-wise dispersion). See
#'   \code{\link[DESeq2]{DESeq}} for more details. (Default="local")
#'
#' @param threads Integer number of threads to run in parallel.
#'
#' @return A list of objects: 
#' 
#'    1st object - Summarized result table merged with anno(PAC) 
#'    
#'    2nd object - Target graphs (p-val distribution and volcano) 
#'    
#'    3rd object - All output from \code{\link[DESeq2]{DESeq}}
#'    
#' @examples
#' 
#' 
#'# Note, these examples will generate some warnings since data is based on
#'# heavily down-sampled fastq files, where many sequences receives low counts
#'# in specific groups.
#'
#'## Load test data
#'load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                 package = "seqpac", mustWork = TRUE))
#'
#'## Simple model with embryonic stages using Wald test with local fit (default)
#'table(pheno(pac)$stage)
#'output_deseq <- suppressWarnings(PAC_deseq(pac, model= ~stage, threads=1))
#'
#'## Batch corrected, graphs are generated for 'stage' (=first in the model)
#'output_deseq <- suppressWarnings(PAC_deseq(pac, model= ~stage + batch,
#'                                           threads=1))
#'
#'## Using pheno_target 
#'output_deseq <- suppressWarnings(PAC_deseq(pac,model= ~stage + batch, 
#'                                           pheno_target=list("batch"),
#'                                           threads=1))
#'
#'## With pheno_target we can change the direction for the comparison
#'# Stage5 vs Stage3 (reverse order):
#'output_deseq <- suppressWarnings(PAC_deseq(pac, model= ~stage + batch, 
#'                          pheno_target = list("stage", c("Stage5", "Stage3")),
#'                          threads=1))  
#'
#'## In the output you find PAC merged results, target plots and output_deseq   
#'names(output_deseq)
#'tibble::as_tibble(output_deseq$result)
#'
#' @importFrom S4Vectors mcols
#' @importFrom stats terms.formula
#' 
#' @export

PAC_deseq <- function(PAC, model, deseq_norm=FALSE, test="Wald", 
                      fitType="local", threads=1, pheno_target=NULL){

## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }    
  cat("\n")
##### Prepare data
  pval <- log2FC <- neglog_padj <- DE <- NULL
  PAC_sub <- PAC
  PAC_sub$Counts <- apply(PAC_sub$Counts, 2, as.integer)
  rownames(PAC_sub$Counts) <- rownames(PAC$Anno)
  
  anno <- PAC_sub$Anno
  pheno <- PAC_sub$Pheno
  
  # Make factors of model columns
  cols <- attr(terms.formula(model), "term.labels")
  cols <- unique(unlist(strsplit(cols, ":")))
  compr<-pheno[,colnames(pheno) %in% cols]
  if(is(compr, "data.frame")){
  if((any(apply(combn(ncol(compr), 2), 2, function(x) identical(compr[, x[1]], compr[, x[2]]))))==TRUE) {
    stop(cat="The column names in model appears to be repeated. \nThis may cause unwanted comparisons. To ensure a correct comparison, please check the colnames in Pheno!")
  }}
  
  for (i in seq.int(length(cols))){
    #Sometime model terms ar complex
    # search for best pheno columns
    trm <- cols[i]
    colnam <- colnames(pheno)
    logi <- NULL
    for(k in seq.int(length(colnam))){
    logi <- c(logi, grepl(colnam[k], trm))
    }
    # Refactor best hit
    colnam_hits <- colnam[logi]
    bst <- nchar(trm) - nchar(colnam_hits)
    log_bst <- bst == min(bst)
    pheno[,colnam_hits[log_bst]]  <- as.factor(pheno[,colnam_hits[log_bst]])
  }
  
  # Prepare pheno target and order factor
  if(is.null(pheno_target)){
      pheno_target <- list(cols[1])
  }
  if(length(pheno_target)==1){
      pheno_target[[2]] <- as.character(unique(PAC$Pheno[,pheno_target[[1]]]))
    }
  
  #Check whether the pheno_target picks up one or many columns (will trg be
  #a factor or a data.frame)  
  trg <- pheno[,colnames(pheno) == pheno_target[[1]]]
  if(is(trg, "data.frame")){
    warning(cat="The values in designated pheno_target may not unique column names. \nThis may cause unwanted comparisons. 
            To ensure a correct comparison, please check the values in Pheno!")
  }
  
  mis <- !levels(trg) %in% pheno_target[[2]] 
  pheno[,colnames(pheno) == pheno_target[[1]]] <- factor(
    trg, levels=c(rev(pheno_target[[2]]),levels(trg)[mis]))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = PAC_sub$Counts, 
                                        colData = droplevels(pheno), 
                                        design=model)
  dds <- DESeq2::estimateSizeFactors(dds)
  
  ### DEseq analysis and extract result table
  BiocParallel::register(BiocParallel::MulticoreParam(workers=threads))
  dds_fit <- DESeq2::DESeq(dds, test=test, fitType=fitType, parallel = TRUE)
  res_nam <- DESeq2::resultsNames(dds_fit)
  if(!is.null(pheno_target)){
     target_nam <- res_nam[grepl(pheno_target[[1]], res_nam)]
     target_nam <- target_nam[1]
  }else{
     target_nam <- res_nam[2]
  }  
  res_DESeq2 <- DESeq2::results(dds_fit, name=target_nam)
  comp <- strsplit(S4Vectors::mcols(res_DESeq2)[2][,1][2], ": ")[[1]][2]
  cat("\n")
  cat("\n")
  cat(paste0("** ", comp, " **"))
  cat("\n")
  
  # Print summary working with different DESeq versions
  test <- try(cat(DESeq2::summary(res_DESeq2)), silent = TRUE)
  if(methods::is(test,"try-error")){ 
    cat(DESeq2::DESeqResults(res_DESeq2))
  }else{
    print(test)
  }
  
  res_DESeq2_df <- as.data.frame(res_DESeq2)
  anno_filt <- anno[match(rownames(res_DESeq2), rownames(anno)),,drop=FALSE]

  if(!identical(rownames(dds), rownames(res_DESeq2))){
      stop("\nNot identical ids in result and dds.",
           "\nHave you been messing with the code?\n")
      }
  if(!identical(rownames(dds), rownames(anno_filt))){
      stop("\nNot identical ids in anno and dds.\n")
      }
  res_counts <- cbind(
    data.frame(feature_ID=rownames(dds)), 
    res_DESeq2_df, 
    DESeq2::counts(dds, normalized=deseq_norm), 
    anno_filt)
  colnames(res_counts)[colnames(res_counts)=="log2FoldChange"] <-  paste(
    "log2FC", res_nam[2], sep="_")
  
  ###  Make plots
  logi_thresh <- ifelse(rowSums(
    cbind(
      res_DESeq2_df$padj <= 0.1,
      res_DESeq2_df$log2FoldChange <=-1 | res_DESeq2_df$log2FoldChange >=1))==2,
    "pass", "not_pass")
    df_plot <- data.frame(
      pval=res_DESeq2_df$pvalue, 
      neglog_padj=-log10(res_DESeq2_df$padj), 
      log2FC=res_DESeq2_df$log2FoldChange, 
      DE=logi_thresh)
  
  p <- ggplot2::ggplot(data=df_plot, ggplot2::aes(x=pval)) + 
    ggplot2::geom_histogram(breaks=seq(0.0, 1.0, by=0.025), 
                            col="black", fill="green", alpha=1) +
    ggplot2::labs(title="p-value distributions", 
                  subtitle=comp, x="p-value", y = "Number of features") +
    ggplot2::theme_classic()
  
  vcano <- ggplot2::ggplot(df_plot, ggplot2::aes(x=log2FC, y=neglog_padj)) +
    ggplot2::geom_hline(yintercept=1, col="black", size=0.1)+
    ggplot2::geom_vline(xintercept=c(-1, 1), col="black", size=0.1)+
    ggplot2::geom_point(ggplot2::aes(colour = DE), size=1) +
    ggplot2::scale_colour_manual(values = c("not_pass"="grey", "pass"= "red")) +
    ggplot2::labs(title="Volcano plot DE features", 
                  subtitle="red points:\nlog2FC >=1  &  p-adj <=0.1",
                  x="Log2 fold changes", y = "-log10 p-value") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0),
                   legend.position = "none")
            
  res_lst <- list(result=res_counts, 
                  plots=list(histogram=p, volcano=vcano), 
                  output_deseq= res_DESeq2)
  print(cowplot::plot_grid(plotlist=res_lst$plots))
  return(res_lst)
}



