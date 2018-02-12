\name{mvGSTsamples}
\docType{data}
\alias{mvGSTsamples}
\alias{obatoclax.pvals}
\alias{parathyroid.pvals}
\alias{obatoclax.mvGST}
\alias{parathyroid.mvGST}
\title{mvGSTsamples:  Sample objects for the mvGST package }
\description{
  This includes the following four objects:
\tabular{ll}{
  \tab \cr
  \code{obatoclax.pvals} \tab a matrix with rows for p-values and columns for contrasts, \cr
  \tab based on GSE36149 data from NCBI GEO \cr
  \code{parathyroid.pvals} \tab a matrix with rows for p-values and columns for contrasts, \cr
  \tab based on parathyroidGenesSE data from \emph{parathyroidSE} package \cr
  \code{obatoclax.mvGST} \tab object returned by \code{profileTable} run on \code{obatoclax.pvals} \cr
  \code{parathyroid.mvGST} \tab object returned by \code{profileTable} run on \code{parathyroid.pvals} \cr
  

  }
}
\usage{data(mvGSTsamples)}
\format{
  This object contains the four objects described above.
}

\value{
  The name of the data set specified.
}

\details{

See \emph{mvGST} package vignette for details on these objects (including how they were constructed): \code{vignette("mvGST")}

The intended use for these objects is to demonstrate the methods coded in the \emph{mvGST} package.
}

\references{ 
Love M., "parathyroidSE: SummarizedExperiment for RNA-Seq of primary culture
parathyroid tumors by Haglund et al., J Clinical Endocrinol Metab 2012." R package
version 1.2.0.

Urtishak K.A., Edwards A.Y., Wang L.S., Hudome A., et al. (2013),
"Potent obatoclax cytotoxicity and activation of triple death mode 
killing across infant acute lymphoblastic leukemia," Blood 121(14):2689-703.

 }

 \examples{
data(mvGSTsamples)
head(obatoclax.pvals)
head(parathyroid.pvals)
obatoclax.mvGST
parathyroid.mvGST

# See package vignette for discussion of these objects: 
#    vignette("mvGST")

}
\keyword{datasets}
