#' Annotate against a GTF file
#'
#' This function will annotate a PAC object using input from a GTF/GFF file.
#'
#' Given a PAC object and a gtf formatted annotation file(s), this function will
#' attempt to annotate sequences mapped to a reference genome against genomic
#' coordinates in the gtf file(s). In case no genomic mapping coordinates are
#' available in the PAC object, the function provides a backdoor into the
#' \emph{PAC reannotation} workflow, where genome mapping is performed using
#' Bowtie.
#' 
#' @family PAC analysis
#'   
#' @seealso \url{https://github.com/OestLab/seqpac} for updates.
#'
#' @param PAC PAC-list object. The Anno object needs to contain genome
#'   coordinates in the following columns: "chr","start","end", and "strand".
#'   These coordinates are eaisly obtained by using \code{PAC_genome}
#'   
#' @param genome Character indicating the path to a reference genome file in
#'   fasta format. Note that a bowtie index, created by for example
#'   Rbowtie::bowtie-build,  must be contained in the same folder and having the
#'   same basename as the fasta.
#'   
#' @param mismatches Integer indicating the number of allowed mismatches for the
#'   genome mapping.
#'   
#' @param return Character indicating what information to return. 
#' 
#'    If return="simplify" (default), a table with one unique
#'    sequence per row will be returned. Multiple hits between genomic
#'    coordinates and gtf coordinates will be merged and only unique annotations
#'    will be reported.
#'    
#'    If return="full", a list will be returned containing 1 table reporting all
#'    annotations for each genomic coordinate of each unique sequence.
#'    
#'    If return="all", both a simplified table and a full annotation list will
#'    be returned as a list.
#'    
#'    If return="merge", a simplified table will be merged with the Anno table
#'    of the provided PAC object, and the updated PAC object containing the new
#'    annotations will be returned.
#'
#' @param gtf Named list of characters, indicating file path(s) to other
#'   gtf files with differing formats. Can also directly be provided as a tibble
#'   dataframe in a named list.
#' 
#' @param targets Named list of character vectors indicating target columns
#'   in the listed gtf files in \emph{gtf_other}. Important, the listed objects
#'   must have the same length and names as in \emph{gtf}. The vector
#'   indicates the column names as if imported by rtracklayer::readGFF.
#'
#' @param stranded Logical whether mapping should be strand specific. If
#'   stranded=TRUE, then hits between feature and PAC read sequence will not be
#'   reported if the feature is located on opposite strand. If
#'   stranded=FALSE (default), then both sense and anti-sense overlaps will be
#'   reported. Note, stranded=FALSE is recommended since PAC_gtf will report on
#'   which strand each feature/sequence is mapping. Thus, strand specific
#'   analysis can be done in hindsight.
#'   
#' @param threads Integer indicating the number of parallel processes.
#'   
#' @return List, tibble dataframe or updated PAC object. See option
#'   \emph{return} for more information.
#'   
#'   
#' @examples
#'
#' # Load PAC
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
#'                  package = "seqpac", mustWork = TRUE))
#' 
#' # Create a gtf file from PAC coordinates
#' anno <- anno(pac)
#' anno <- anno[!grepl("Warning", anno$mis0_chromosomes_genome),]
#' anno <- anno[!is.na(anno$mis0_chromosomes_genome),]
#' coord <- anno$mis0_chromosomes_genome
#' coord <- suppressWarnings(do.call("rbind", strsplit(coord, "\\|"))[,1])
#' coord <- suppressWarnings(do.call("rbind", strsplit(coord, "\\;start=|;")))
#' gr <- GenomicRanges::GRanges(seqnames=coord[,1], 
#'                              IRanges::IRanges(as.numeric(coord[,2]), 
#'                                               as.numeric(coord[,2])+
#'                                               anno$Size ), 
#'                              strand=coord[,3])
#' 
#' GenomicRanges::mcols(gr) <- data.frame(biotype=anno$Biotypes_mis3, 
#'                                        bio_zero=as.character(anno$mis0_bio))
#' spl <- sample(1:length(gr), round(length(gr)*0.3), replace=FALSE)
#' gr1 <- gr[spl]
#' gr2 <- gr[!1:length(gr) %in% spl]
#' 
#' # Prepare temp folder and save artifical gtf
#' out1 <- paste0(tempdir(), "/temp1.gtf")
#' out2 <- paste0(tempdir(), "/temp2.gtf")
#' 
#' rtracklayer::export(gr1, out1, format="gtf")
#' rtracklayer::export(gr2, out2, format="gtf")
#' 
#' # Make sure PAC contains full genome coordinates
#' # (In the add_reanno function you may choose how many coordinates to report)
#' # (If there are more, a 'Warning' will be added to the annotation)
#' # (Here we remove those to avoid problems)
#' 
#' new_anno <-  anno(pac) [, grepl("chromosomes_genome|Size", 
#'                                colnames(anno(pac)))]
#' test <- new_anno[, 3:ncol(new_anno)]
#' test <- apply(test, 2, function(x){substr(x, 1, 7)})
#' test <- apply(test, 1, function(x){paste(x, collapse = "")})
#' new_anno$temp <- ifelse(grepl("Warning",  test), "rm", "keep")
#' anno(pac) <- new_anno
#' pac_small <- PAC_filter(pac, subset_only=TRUE, size = c(20,22),
#'                   anno_target=list("temp", "keep"))
#'
#' #  Run PAC_gtf
#' gtf <- list(gtf1=out1, gtf2=out2)
#' target <- list(gtf1=c("biotype","bio_zero"), gtf2=c("biotype","bio_zero"))
#' 
#' pac_merge <- PAC_gtf(pac_small, mismatches=0, return="merge", 
#'                     gtf=gtf, target=target, threads=2)
#' 
#' 
#' @export


PAC_gtf<- function(PAC, genome=NULL, mismatches=3, return="simplify", 
                   stranded=FALSE, gtf=NULL, targets=NULL, 
                   threads=1){
  ## Check S4
  if(isS4(PAC)){
    tp <- "S4"
    PAC <- as(PAC, "list")
  }else{
    tp <- "S3"
  }
  ##### Setup general ####################################
  seqs <- rownames(PAC$Anno)
  
  ##### Import gtfs ####################################
  if(!is.list(gtf)){
    gtf <- list(gtf)
    if(is.null(names(gtf))){
      names(gtf) <- seq.int(length(gtf))
    }
  }
  cat("\nImporting gtf files ...")
  gtf_lst <- lapply(gtf, function(x){
    if(is.character(x[1])){    
      if(file.exists(x)){
        gtf <- tibble::as_tibble(rtracklayer::readGFF(x), 
                                 .name_repair="minimal")
        return(gtf)
      }else{
        return("No gtf-formatted file was found.")  
      }
    }
    return(x)
  })
  
  
  ##### Check gtf file and columns ####################################  
  for (i in seq.int(length(gtf_lst))){
    nam <- names(gtf_lst)[i]
    logi_fl <-  methods::is(gtf_lst[[i]], "tbl_df")
    if(!logi_fl){
      stop("
           \nFile was missing for '", nam, "' gtf. Please, check file path.")
    }
    
    logi_col <- sum(
       names(gtf_lst[[i]]) %in% c("seqid", "start", "end", "strand"))== 4
    if(!logi_col){
      stop("
           \nInput gtf '", nam, "' does not contain essential columns",
           "\n('seqid', 'start', 'end', 'strand'). Please, check column",
           "\nnames using for example rtracklayer::readGFF('<path_to_gtf>').")
    }
    logi_target <- sum(
       names(gtf_lst[[i]]) %in% targets[[i]])== length(targets[[i]])
    miss_trg <- which(!targets[[i]] %in% names(gtf_lst[[i]]))
    miss_trg <- paste(targets[[i]][miss_trg], collapse="; ") 
    if(!logi_target){
      stop("
           \nInput gtf '", nam, "' does not contain all target columns.",
           "\nPlease, check correct names for column(s) ", miss_trg,
           "\nusing for example rtracklayer::readGFF('<path_to_gtf>').")
    }
  }
  ##### Setup genome and run reanno if necessary ########################
  # If user do not know columns
  if(is.null(genome)){
    anno_genome <- tibble::as_tibble(PAC$Anno, .name_repair="minimal")
  }else{
    # If user know the columns:
    if(any(!file.exists(genome))){
      cat("\nNo genome fasta was found, will atempt", 
          "\nfinding genome columns in PAC@Anno.")
      anno_genome <- tibble::as_tibble(PAC$Anno[, genome, drop=FALSE], 
                                       .name_repair="minimal")
    }else{
      # If user provide fasta reference:
      if(file.exists(genome)){
        cat("\nInput genome was an existing file. Will treat it as a",
            "\nfasta reference and make a denovo reannotation using bowtie.",
            "\nSee ?map_reanno or ?vingette for details.\n")
        outpath <- tempfile(pattern = "", fileext = "")
        # Convert to windows format
        if(grepl("win|WIN|Win", Sys.info()["sysname"])){
            outpath <- gsub("\\", "/", outpath, fixed=TRUE)
        }
        err <- try(map_reanno(PAC, ref_paths=list(genome=genome),
                                 output_path=outpath, type="internal",
                                 mismatches=mismatches, override=TRUE,
                                 import ="genome", threads=threads),
                   silent = TRUE)
        if(!is.null(err)){
        outpath <- tempfile(pattern = "", fileext = ".out")
        err2 <- try(map_reanno(PAC, ref_paths=list(genome=genome), 
                              output_path=outpath, type="external", 
                              mismatches=mismatches, override=TRUE,
                              import ="genome", threads=threads), silent = TRUE)
        if(!is.null(err2)){
            stop(
                "\nFunction map_reanno failed. Possible reasons:",
                "\n\tNo bowtie installation",
                "\n\tBad fasta reference",
                "\n\tNo bowtie index (see ?Rbowtie::bowtie_build)",
                "\n\nLast log says:\n", err2)
          }
        }
        reanno <- make_reanno(outpath, PAC=PAC, 
                              mis_fasta_check=TRUE, output="list" )
        anno_genome <- add_reanno(reanno, type="genome", 
                                  mismatches = mismatches, genome_max="all")
        rm(reanno)
        unlink(outpath, recursive = TRUE)
      }
    }
  }
  # Extract genome columns
  srch <- paste0("^mis", 0:mismatches)
  c_nams <- colnames(anno_genome)
  logi_col1 <- grepl("_genome$", c_nams)
  logi_col2 <- grepl(paste(srch, collapse="|"), c_nams)
  nam_col <- gsub("mis0_", "", c_nams[grepl("^mis0_", c_nams)])
  mis_col <- grepl(paste0("^", nam_col), c_nams)
  logi_coord <- rowSums(cbind(logi_col1, logi_col2))==2
  try_err <-  try(prefix <- do.call("rbind", 
                                   strsplit(c_nams[logi_coord], "_"))[,1], 
      silent=TRUE)
  if(!is(try_err, "try-error")){
      warning("Failed to generate prefix")
      }
  
  if(any(duplicated(prefix))){
    stop(
      "\nFound more than 1 possible genome annotation.",
      "\nAutomatic column identification only works when the results",
      "\nof one fasta reference genome is reported in PAC@Anno.",
      "\nPlease, specify in 'genome=' the exact columns reporting the",
      "\nalignments against 1 reference genome (column prefix: mis0_, mis1_,",
      "\nmis2_ etc.) or proved a path to a bowtie indexed fasta genome file.")
  }
  coord_genome <- anno_genome[,logi_coord, drop=FALSE]
  mis_genome <- anno_genome[,mis_col, drop=FALSE]
  # Check that all columns are included
  logi_check  <- apply(coord_genome, 2, function(x){
    any(grepl("Warning>\\d", x))
  })
  if(any(logi_check)){
    stop(
    "\nCoordinates in genome appear to have truncation marks ('[W]arning>').",
      "\nThis indicates that genome was mapped reporting a limited number",
      "\nof alignments per sequence. Please, rerun the 'add_reanno' function",
      "\nusing max_genome='all', or provide the path to a fasta genome",
      "\nreference in 'genome'.")
  }
  
  ##### Organize mapped coordinates ####################################
  cat("\n\nReorganizing coordinates ...")
  coord_lst <- list(NULL)
  lng_all <- PAC$Anno[,colnames(PAC$Anno) %in% c("Length", "Size"), drop=TRUE]
  mis_incl <- unique(mis_genome[[1]])
  mis_incl <- mis_incl[mis_incl %in% paste0("mis", 0:mismatches)]
  for(i in seq.int(nrow(mis_genome))){
    mis <- mis_genome[[1]][i]
    lng <- lng_all[i]
    if(!mis %in% mis_incl){
      tb <- tibble::tibble(seqid=NA, start=NA, end=NA, strand=NA)
    }else{  
      targt <- coord_genome[i, grepl(mis, names(coord_genome))]
      splt <- unlist(strsplit(targt[[1]], "\\|"))
      tb <- tibble::as_tibble(do.call("rbind", strsplit(splt, ";")), 
                              .name_repair="minimal")
      tb[[2]] <- as.numeric(gsub("start=", "", tb[[2]]))
      tb <- tibble::tibble(seqid=tb[[1]], start=tb[[2]], 
                           end=tb[[2]]+lng, strand=tb[[3]])
    }
    coord_lst[[i]] <- tb
    names(coord_lst)[i] <- mis
  }
  
 
  ##### Check all chromosomes names in gtf ####################################
  cat("\n\nChromosome names compatibility check ...")
  chrm_nam <- unique(unlist(lapply(coord_lst, function(x){x$seqid})))
  chrm_nam <- chrm_nam[!is.na(chrm_nam)]
  tot_chrom <- length(chrm_nam)
  
  for(i in seq.int(length(gtf_lst))){
    chrom_in_gtf <- sum(chrm_nam %in% unique(gtf_lst[[i]]$seqid))
    if(chrom_in_gtf==0){
      stop("\nYour chromosome names in 'gtf_", 
           names(gtf_lst)[i], 
           "' input \ndid not match your genome alignments. Did you forget",
           "\nto convert between database formats (e.g. UCSC vs Ensembl vs",
           "\nNCBI)? Please, check ?PAC_gtf and vingette(sepac) to", 
           "\ninstructions on how to harmonize your gtf files.")
    }
    if(chrom_in_gtf/tot_chrom<0.05 ){
      cat("\n")
      warning(
        "Low overlap between chromosome names of 'gtf_", 
        names(gtf_lst)[i], 
        "'\ninput and genome alignments was detected. Please, double",
        "\ncheck that conversion between database formats (e.g. UCSC vs",
        "\nEnsembl vs NCBI) was completed. Check ?PAC_gtf and vingette(sepac)",
        "\nfor instructions on how to harmonize your gtf files to",
        "\nthe reference genome (Note - this may also depend on short gtf file).")
    }
    if(!chrom_in_gtf/tot_chrom==1 ){
      cat("\n")
      cat(
        paste0(
          "  -- Note, not all chrosomsome names were represented in 'gtf_", 
          names(gtf_lst)[i],
          "'.\n     Could be wise to double check that chrosomsome names ",
          "\n     are compatible."))
    }
  }
  
  ##### Annotating ####################################
  cat("\n\nAnnotating against the gtf file(s)  ...")
  # Convert gtf to gr
  all_chrom_nam <- unique(c(as.character(
    unique(unlist(lapply(gtf_lst, function(x){
      x$seqid
      })))), chrm_nam))
  gtf_gr<- lapply(gtf_lst, function(x){
    gr <- GenomicRanges::GRanges(seqnames=as.character(x$seqid), 
                                 IRanges::IRanges(start=as.integer(x$start), 
                                                  end=as.integer(x$end)), 
                                 strand=as.character(x$strand))
    GenomeInfoDb::seqlevels(gr) <- all_chrom_nam
    return(gr)
  })
  
  # Convert genome coordinates to gr
  `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(threads) # Do not use parallel::makeClusters!!!
  coord_gr   <- foreach::foreach(i=seq.int(length(coord_lst)), 
                                 .packages=c("GenomicRanges"), 
                                 .final = function(x){
                                   names(x) <- names(coord_lst); return(x)
                                   }) %dopar% {
    x <- coord_lst[[i]]
    if(is.na(x[[1,1]])){
      gr <- NA
    }else{
      if(stranded=="TRUE"){
        gr <- GenomicRanges::GRanges(seqnames=x$seqid, 
                                     IRanges::IRanges(start=x$start,
                                                      end=x$end), 
                                     strand=x$strand)
        GenomeInfoDb::seqlevels(gr) <- all_chrom_nam
      }
      if(stranded=="FALSE"){
        gr <- GenomicRanges::GRanges(seqnames=x$seqid, 
                                     IRanges::IRanges(start=x$start, 
                                                      end=x$end), 
                                     strand=NA)
        GenomeInfoDb::seqlevels(gr) <- all_chrom_nam
      }
    }
    return(gr)
  }
  doParallel::stopImplicitCluster()
  
  # Extraction loop
  doParallel::registerDoParallel(threads)
  anno_lst <- list(NULL)
  for(i in seq.int(length(gtf_gr))){
    gtf_nam <- names(gtf_gr)[i]
    cat(paste0("\n   |--> Extract and compile '", gtf_nam, "' ..."))  
    trg_cols <- targets[[i]]

    # Run overlap and extract anno 
    coord_anno <-  foreach::foreach(t=seq.int(length(coord_gr)), 
                                    .packages=c("GenomicRanges"), 
                                    .final = function(t){
                                      names(t) <- names(coord_gr); return(t)
                                      }) %dopar% {
      x <- coord_gr[[t]]
      if(paste(x[1])=="NA"){
        uni_gtf <- tibble::as_tibble(matrix(NA, nrow=1, 
                                            ncol=length(trg_cols)), 
                                     .name_repair="minimal")
        names(uni_gtf) <- trg_cols
      }else{
        try_err<- try(olap <- GenomicRanges::findOverlaps(x, gtf_gr[[i]]),
                      silent=TRUE)
        if(is(try_err, "try-error")){
          warning("Something went wrong in GenomicRanges::findOverlaps")
        }
                gtf <- gtf_lst[[i]][S4Vectors::subjectHits(olap), trg_cols]
        
        uni_gtf <- data.frame(matrix(NA, nrow=length(x), 
                                     ncol=length(trg_cols)))
        names(uni_gtf) <- trg_cols
        if(!nrow(gtf)==0){
          splt <- split(gtf, as.factor(S4Vectors::queryHits(olap)))
          uni_anno <- lapply(splt, function(z){
            uni <- apply(z, 2, function(y){
              y <- y[!is.na(y)]
              paste(sort(unique(y)), collapse="|")
            })
            return(uni)
          })
          uni_anno <- do.call("rbind", uni_anno)
          uni_gtf[rownames(uni_anno),] <- uni_anno
        }
      }
      return(uni_gtf)
    }
    anno_lst[[i]] <- coord_anno
    names(anno_lst)[i] <- gtf_nam
  }
  doParallel::stopImplicitCluster() 
  
  ##### Returning objects ####################################  
  if(return %in% c("full","all")){
    cat("\n\n")
    full_lst <- list(NULL)
    loop <- length(coord_lst)
    for(i in seq.int(loop)){
      cat("\rGenerating full annotation list: ", i, "/", loop)
      utils::flush.console() 
      tab <- coord_lst[[i]]
      for(z in seq.int(length(anno_lst))){
        tab <- cbind(tab, anno_lst[[z]][[i]])
      }
      full_lst[[i]] <- tibble::as_tibble(tab, .name_repair="minimal")
      names(full_lst)[i] <- paste(seqs[i], names(coord_lst)[i], sep="_")
    }
  }
  
  if(return %in% c("simplify", "all", "merge")){
    cat("\n\n")
    simp_lst <- list(NULL)
    loop <- nrow(coord_genome)
    for(i in seq.int(loop)){
      cat("\rGenerating simplified annotation table: ", i, "/", loop)
      utils::flush.console() 
      tab <- coord_genome[i,]
      for(z in seq.int(length(anno_lst))){
        simp_tab <- apply(anno_lst[[z]][[i]], 2, function(x){
          simp <- sort(unique(unlist(strsplit(as.character(x), "\\|"), 
                                     use.names = FALSE)))
          simp <- paste(simp[!is.na(simp)], collapse="|")
          if(nchar(simp) == 0){
            simp <- NA
          }
          return(simp)
        })
        tab <- cbind(tab, t(as.data.frame(simp_tab)))
      } 
      simp_lst[[i]] <- tab
      names(simp_lst)[i] <- seqs[i]
    }
    simp_df <- tibble::as_tibble(cbind(data.frame(seqs=as.character(seqs)),
                                       mis_genome, do.call("rbind", simp_lst)), 
                                 .name_repair="minimal")
  }
  
  if(return=="simplify"){
    return(simp_df)
  }
  if(return=="full"){
    return(full_lst)
  }
  if(return=="all"){
    return(list(simplify=simp_df, full=full_lst))
  }
  if(return=="merge"){
    stopifnot(identical(rownames(PAC$Anno), as.character(simp_df$seqs)))
    PAC$Anno <- cbind(PAC$Anno, simp_df[,-1])
    stopifnot(PAC_check(PAC))
    if(tp=="S4"){
      return(as.PAC(PAC))
    }else{
      return(PAC)
    }    
  }
}
