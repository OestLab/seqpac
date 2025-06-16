#' Wrapper for annotation of PAC object
#'
#' \code{PAC_annotate} Annotates the provided PAC object to a genome or biotypes
#'
#' Given the path to either a reference genome, a biotype reference or both, PAC_annotate will wrap 
#' around all reanno-functions (\code{\link{map_reanno}}, \code{\link{make_reanno}}, and
#' \code{\link{add_reanno}}) to annotate the PAC object. This wrapper comes with the caveat that
#' it uses default values in most functions and is finetuned towards basic biotype annotation
#' against the Ensembl ncRNA fasta reference. If that does not suit your analysis, we recommend
#' you to run each function on its own for its full potential!
#'
#' @family PAC reannotation
#'
#' @seealso \url{https://github.com/OestLab/seqpac} for updates on the current
#'   package.
#'
#' @param genome Character indicating path to reference genome in fasta (.fa) format to use
#' for annotation by \code{\link{map_reanno}} with import="genome". Reference genome should have a
#' bowtie index, please see ??Rbowtie::bowtie_build for instructions. Default=NULL, where no genome-based
#' mapping will be performed.
#'
#' @param biotype Character indicating path to reference fasta (.fa) file for bioinformatic annotation
#' by \code{\link{map_reanno}} with import="biotype". This wrapper is justerad to the ncRNA fasta from Ensembl,
#' with biotypes such as rRNA, tRNA and miRNA available. Reference fasta should have a
#' bowtie index, please see ??Rbowtie::bowtie_build for instructions. Within this wrapper, a hierarcy is made with
#' \code{\link{simplify_reanno}}, with the following hiearchy: rRNA, miRNA, tRNA, snoRNA, snRNA. If other references wish
#' to be used, please run each function separately. Default=NULL, where no biotype-based
#' mapping will be performed.
#'   
#' @param pac PAC-object used for mapping.
#'    
#' @return PAC object with annotation information retrieved from mapping
#'   
#' @examples
#' 
#' ###########################################################
#'
#' genome <- system.file("extdata/mycoplasma_genome", "mycoplasma.fa",
#'                      package = "seqpac", mustWork = TRUE)
#' genome_dir <- dirname(genome)
#' 
#' if(!sum(stringr::str_count(list.files(genome_dir), ".ebwt")) ==6){
#'   Rbowtie::bowtie_build(genome,
#'                         outdir=genome_dir,
#'                         prefix="mycoplasma", force=TRUE)
#' }
#'
#' ##  load a PAC-object and remove previous mapping from anno:
#' load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata",
#'                  package = "seqpac", mustWork = TRUE))
#' anno(pac) <- anno(pac)[,1, drop = FALSE]
#' 
#' #as we only test this against the mycoplasma genome,
#' #we expect a small overlap (0.03%)
#' 
#' pac <- PAC_annotate(pac=pac,
#'                     genome = genome)
#'
#' @export



PAC_annotate <- function(genome=NULL,
                         biotype=NULL,
                         pac){
  
  cat("Using the working directory as output folder: ",getwd())
  
  # if(is.null(reference)){
  #   cat("No reference fasta were given to the function. Will now download ")
  # }
  # 
  
  if(!is.null(genome)){
    
    #    here we could make a function that picks up a reference online ... 
    # if(!file.exists(genome)){
    #   cat("No reference fasta were given to the function. Will now download reference genome for
    #       mmu () to use instead... ")
    #   
    #   genome <- download.file("https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz")
    # }
    #if(){} - here, check whether a file or a website to use like biomart or bsgenome
    #perform mapping with "genome" functionailty
    
    map_reanno(pac, import="genome", ref_paths=list(genome=genome), 
               output_path=getwd(), mismatches=0, override = TRUE)
    reanno<-make_reanno(getwd(), PAC=pac, mis_fasta_check=TRUE)
    pac <- add_reanno(reanno, type="genome", mismatches=0, merge_pac=pac)
    
  }
  
  if(!is.null(biotype)){
    #perform mapping with "biotype" functionality
    map_reanno(pac, import="biotype", ref_paths=list(biotype=biotype), 
               output_path=getwd(), mismatches=0, override=TRUE)
    reanno<-make_reanno(getwd(), PAC=pac, mis_fasta_check=TRUE)
    bio_search <- list(biotype=c("rRNA", "tRNA", "miRNA",
                                 "snoRNA", "snRNA", "piRNA"))
    pac <- add_reanno(reanno, type="biotype", mismatches=0, merge_pac=pac, bio_search=bio_search,
                      bio_perfect=FALSE)
    pac <- simplify_reanno(input=pac, mismatches=0,
                           merge_pac=TRUE,
                           hierarchy=list(rRNA="biotype_rRNA",
                                          miRNA="biotype_miRNA",
                                          tRNA="biotype_tRNA",
                                          snoRNA="biotype_snoRNA",
                                          snRNA="biotype_snRNA"))
    
  }
  
  return(pac)
  
}
