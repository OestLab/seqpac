% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAC_summary.R
\name{PAC_summary}
\alias{PAC_summary}
\title{Summarizes PAC objects}
\usage{
PAC_summary(
  PAC,
  norm = "counts",
  type = "means",
  pheno_target = NULL,
  rev = FALSE,
  merge_pac = TRUE
)
}
\arguments{
\item{PAC}{PAC object containing a Pheno data.frame with samples as row
names and a Counts table with raw counts. Optionally, the PAC object may
contain normalized counts tables, saved in the norm list ('folder'). Such
normalized table can be generated using the \code{\link{PAC_norm}}
function.}

\item{norm}{Character indicating what type of data to be used. If 'counts',
the raw counts in Counts will be used (default). Given any other value, the
function will search for the value as a name on a data.frame stored in the
normalized list-folder.}

\item{type}{Character indicating what type of summary to be applied to the
data. The function currently supports:
type="means"         # Group means
type="sd"            # Group standard deviation
type="se"            # Group standard error of the mean
type="log2FC"        # Group log2 fold changes against other groups
type="log2FCgrand"   # Group log2 fold changes against a grand mean 
type="percentgrand"  # Group log2 fold changes against a grand mean}

\item{pheno_target}{List with: 1st object being a character vector
of target column in Pheno 2nd object being a character vector of the target
group(s) in the target Pheno column (1st object).}

\item{rev}{Logical whether pairwise comparisons (e.g. log2FC) should be
reversed (default=FALSE).}

\item{merge_pac}{Logical whether simplified annotation vector should
automatically be added to the Anno object of the input PAC list object
(default=TRUE). If \code{merge_pac=FALSE} a data.frame is returned.}
}
\value{
A PAC object with a pheno_summary folder containing the summarized
  data in a data.frame. The data.frame will be named according to the
  pheno_target, type and norm input.
}
\description{
\code{PAC_summary} summarizes data stored in a PAC object.
}
\details{
Given a PAC object this function summarize data in counts(PAC) or in the norm
'folder' according to a grouping column in pheno(PAC).
}
\examples{

load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                 package = "seqpac", mustWork = TRUE))

PAC_check(pac) # TRUE

# Easy to generate simple group summaries 
pac <- PAC_summary(pac, norm = "cpm", 
                   type = "means", pheno_target=list("stage"))       
pac <- PAC_summary(pac, norm = "cpm", 
                   type = "se", pheno_target=list("stage"))
pac <- PAC_summary(pac, norm = "cpm", 
                   type = "log2FC", pheno_target=list("stage"))

names(summary(pac))               # Names of individual summaries
head(summary(pac)$cpmMeans_stage) # View individual summaries
summary(pac)  # View merge summaries
df <- as.data.frame(tibble::as_tibble(summary(pac))) # Merge multiple summaries
head(df)


# If a pheno_target is left out, a mean of all samples will be returned:
load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                 package = "seqpac", mustWork = TRUE))
pac <- PAC_summary(pac, norm = "cpm", type = "mean")  
pac <- PAC_summary(pac, norm = "cpm", type = "percentgrand")
names(summary(pac))
summary(pac)   

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
\code{\link{PAC_nbias}()},
\code{\link{PAC_norm}()},
\code{\link{PAC_pca}()},
\code{\link{PAC_pie}()},
\code{\link{PAC_saturation}()},
\code{\link{PAC_sizedist}()},
\code{\link{PAC_stackbar}()},
\code{\link{PAC_trna}()},
\code{\link{as.PAC}()},
\code{\link{filtsep_bin}()},
\code{\link{map_rangetype}()},
\code{\link{tRNA_class}()}
}
\concept{PAC analysis}
