% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PAC_check.R
\name{PAC_check}
\alias{PAC_check}
\title{Checks PAC format}
\usage{
PAC_check(PAC, sample_check = NULL, seq_check = NULL)
}
\arguments{
\item{PAC}{PAC-list object.}

\item{sample_check}{Can be provided either as a data frame or the name of a
data frame in PAC. Will check both the row and column names of the object
for identical sample names in PAC. A simple character vector can also be
provided.}

\item{seq_check}{Can be provided either as a data frame or the name of an
data frame in PAC. Will check only the row names of the object for
identical sequence names in PAC. A simple character vector can also be
provided.}
}
\value{
If PAC format is correct will return logical TRUE, else it will throw
  an error.
}
\description{
\code{PAC_check} Quickly checks rownames in PAC.
}
\details{
Given a PAC object the function will cross check the rownames and colnames
across Pheno, Anno and Counts.
}
\examples{

load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", 
                  package = "seqpac", mustWork = TRUE))

PAC_check(pac)  # TRUE

# Remove first row in Anno and check compatibility
new_Anno <- anno(pac)[-1,]
#PAC_check(pac, seq_check=new_Anno) # Error

# Add to pac an check
# anno(pac) <- new_Anno #error

}
\seealso{
\url{https://github.com/OestLab/seqpac} for updates on the current
  package.

Other PAC generation: 
\code{\link{PAC_create}()},
\code{\link{make_PAC}()},
\code{\link{make_counts}()},
\code{\link{make_cutadapt}()},
\code{\link{make_pheno}()},
\code{\link{make_trim}()},
\code{\link{merge_lanes}()}
}
\concept{PAC generation}
