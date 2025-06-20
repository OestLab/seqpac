% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_reanno.R
\name{import_reanno}
\alias{import_reanno}
\title{Imports annotation from reannotation mapping}
\usage{
import_reanno(
  bowtie_path,
  threads = 1,
  coord = FALSE,
  report = "minimum",
  reduce = NULL
)
}
\arguments{
\item{bowtie_path}{Path to a directory where bowtie output files can be
found.}

\item{threads}{Integer stating the number of parallel jobs.}

\item{coord}{Logical whether or not mapping coordinates should be reported
when report="full".}

\item{report}{Character vector indicating what to report "minimum" or "full"
(default="minimum").}

\item{reduce}{Character indicating a reference name (without file type
extension) that should be exempted from report="full" and instead will
generate a minimum report.}
}
\value{
List of data frames with additional information from reannotation
  files generated with bowtie. If \emph{report="minimum"}, the function will
  report each hit only as the number of mismatches for a given reference
  file. If \emph{report="full"} the full name reported in the fasta file used
  as reference in the bowtie reannotation will be reported. If a reference
  name is specified in \emph{reduce}, this reference is excerpted from
  \emph{report="full"} and is instead reported as \emph{report="minimum"}.

  Caution: Large references with lots of redundancy (such as pirBase in some
  species) will cause massive character strings if \emph{report="full"} with
  no restrictions. Specifying such references in
  \emph{reduce=<reference_names>} will circumvent this problem.
}
\description{
This function imports bowtie output and summarizes the content.
}
\details{
Given the path to alignment outputs from bowtie, \code{\link{import_reanno}}
will attempt to read these files into R and generate a list of unsorted
data.frames where each row summarizes all annotations for a given sequence.
It is called by \code{\link{map_reanno}} to generate summarized data from
bowtie output saved as .Rdata files. These files can then be re-imported into
R where information is organized using \code{\link{make_reanno}} and finally
added to the annotation table (PAC$Anno) of the original PAC-list object
using \code{\link{add_reanno}}.
}
\section{Important}{
 Re-annotation must have been done using Bowtie default
  output settings, where output files specifically have been named '.txt'.
  The function will not work with other formats such as SAM/BAM formats.

 The basenames of the bowtie output will also be used to
  annotate. Example: bowtie -v 3 -a -f
  '<piRBase.fasta>' '<master_anno.fasta>' piRNA.txt Will annotate sequences
  as piRNA since the txt-file output was named 'piRNA'
}

\examples{

######################################################### 
##### Test import_reanno

### First, if you haven't already generated Bowtie indexes for the included
# fasta reference you need to do so. If you are having problem see the small
# RNA guide (vignette) for more info.

## tRNA reference:
trna_file <- system.file("extdata/trna", "tRNA.fa", 
                         package = "seqpac", mustWork = TRUE)
trna_dir<- gsub("tRNA.fa", "", trna_file)

if(!sum(stringr::str_count(list.files(trna_dir), ".ebwt")) ==6){
  Rbowtie::bowtie_build(trna_file, 
                        outdir=trna_dir, 
                        prefix="tRNA", force=TRUE)
}
##  Then load a PAC-object and remove previous mapping from anno:
load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                 package = "seqpac", mustWork = TRUE))

ref_paths <- list(trna= trna_file)

##  You may add an output path of your choice, but here we use a temp folder:
output <- paste0(tempdir(),"/seqpac/test")

##  Then map the PAC-object against the fasta references. Warning: if you use
# your own data, you may want to use override=FALSE, to avoid deleting previous
# mapping by mistake. keep_temp=TRUE can be used to run import_reanno
# independently.

map_reanno(pac, ref_paths=ref_paths, output_path=output,
           type="internal", mismatches=2,  import="biotype", 
           threads=2, keep_temp=TRUE, override=TRUE)

reanno1 <- import_reanno(output, report="Biotype",  threads=1)

}
\seealso{
\url{http://bowtie-bio.sourceforge.net/index.shtml} for information
  about Bowtie and for Rbowtie:
  \url{https://www.bioconductor.org/packages/release/bioc/html/Rbowtie.html}.
  \url{https://github.com/OestLab/seqpac} for updates on the current package.

Other PAC reannotation: 
\code{\link{PAC_annotate}()},
\code{\link{add_reanno}()},
\code{\link{as.reanno}()},
\code{\link{make_conv}()},
\code{\link{make_reanno}()},
\code{\link{map_reanno}()},
\code{\link{simplify_reanno}()}
}
\concept{PAC reannotation}
