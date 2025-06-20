% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAC_nbias.R
\name{PAC_nbias}
\alias{PAC_nbias}
\title{Generates a nucleotide bias analysis from a PAC object}
\usage{
PAC_nbias(
  PAC,
  position = 1,
  norm = NULL,
  range = NULL,
  anno_target = NULL,
  pheno_target = NULL,
  summary_target = NULL,
  colors = NULL,
  ymax = NULL,
  data_only = FALSE
)
}
\arguments{
\item{PAC}{PAC-list object containing an Anno data.frame with sequences as
row names and a count table with raw counts.}

\item{position}{Integer indicating the nucleotide position from 3' to 5'
position (default=1).}

\item{norm}{Character indicating what type of data to be used. If
type="counts" the plots will be based on the raw Counts. If type="cpm" the
analysis will be done on cpm values returned from \code{PAC_norm} function
and stored in the norm folder of the PAC-list object. The name of any other
table in the norm(PAC) folder can also be used.}

\item{range}{Integer vector indicating the sequence size range
(default=c(min, max)).}

\item{anno_target}{List with: 
1st object being character vector of target
column(s) in Anno, 2nd object being a character
vector of the target biotypes(s) in the target
column (1st object). (default=NULL)}

\item{pheno_target}{List with: 
1st object being character vector of target
column(s) in Pheno, 2nd object being a character
vector of the target group(s) in the target column
(1st object). (default=NULL)}

\item{summary_target}{List with: 
1st object being character vector of target object
in summary(PAC), 2nd object being a character vector
of the target column(s) in the target summary object
(1st object). (default=NULL)}

\item{colors}{Character vector with RGB color codes to be parsed to ggplot2.}

\item{ymax}{Integer that sets the maximum y for all plots (all plots gets the
same y-axes). If ymax=NULL, then ggplot2 will automatically set ymax for
each plot individually (different y-axes).}

\item{data_only}{logical. If data_only=TRUE a data.frame a simple Anno object
is returned with a Size and a Nucleotide bias column. As default,
data_only=FALSE then graphs are returned in addition to data.}
}
\value{
A list of objects: 
              1st object (Histograms::Samples): Individual histograms showing
              the nucleotide ratios per sample over the specified range. 2nd
              object (Data::Samples): Data used to generate the plots.
}
\description{
\code{PAC_nbias} analyses nucleotide bias.
}
\details{
Given a PAC object the function will attempt to extract the ratios of
specific nucleotides at a given position in sequences in the Anno data.frame
in relation to the sequence counts in Counts.
}
\examples{


# Load a PAC-object 
load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                  package = "seqpac", mustWork = TRUE))
output_nbias <- PAC_nbias(pac)
cowplot::plot_grid(plotlist=output_nbias$Histograms)

# Only miRNA (Oops, heavy T-bias on 1st nt; are they piRNA?)  
table(anno(pac)$Biotypes_mis0)
output_nbias <- PAC_nbias(pac, anno_target = list("Biotypes_mis0", "miRNA") )
cowplot::plot_grid(plotlist=output_nbias$Histograms)

# Switch to 10:th nt bias 
output_nbias <- PAC_nbias(pac, position=10, 
                          anno_target = list("Biotypes_mis0", "miRNA"))
cowplot::plot_grid(plotlist=output_nbias$Histograms)

# Summarized over group cpm means
pac_test <- PAC_summary(pac, norm = "cpm", type = "means", 
                        pheno_target=list("stage"), merge_pac=TRUE)
output_nbias <- PAC_nbias(pac_test, summary_target = list("cpmMeans_stage") )
cowplot::plot_grid(plotlist=output_nbias$Histograms)


}
\seealso{
\url{https://github.com/OestLab/seqpac} for updates on the current
  package.

Other PAC analysis: 
\code{\link{PAC_covplot}()},
\code{\link{PAC_deseq}()},
\code{\link{PAC_filter}()},
\code{\link{PAC_filtsep}()},
\code{\link{PAC_gtf}()},
\code{\link{PAC_jitter}()},
\code{\link{PAC_mapper}()},
\code{\link{PAC_norm}()},
\code{\link{PAC_pca}()},
\code{\link{PAC_pie}()},
\code{\link{PAC_saturation}()},
\code{\link{PAC_sizedist}()},
\code{\link{PAC_stackbar}()},
\code{\link{PAC_summary}()},
\code{\link{PAC_trna}()},
\code{\link{as.PAC}()},
\code{\link{filtsep_bin}()},
\code{\link{map_rangetype}()},
\code{\link{tRNA_class}()}
}
\concept{PAC analysis}
