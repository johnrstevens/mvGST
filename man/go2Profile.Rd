\name{go2Profile}
\alias{go2Profile}
\title{
Creates tables showing profiles of specific gene sets.
}
\description{
Creates tables similar to the results tables of an \code{mvGST} object except for one gene set at a time. For each gene set selected, a table is produced with a single one and all zeroes in each column.
}
\usage{
go2Profile(names, object)
}
\arguments{
  \item{names}{
A character vector with the names, ID's, of the gene sets of interest. If the gene set names were not provided by the user, then this should be the GO ID's of the gene sets of interest.
}
  \item{object}{
A \code{mvGST} object with a final \code{results.table}.
}
}

\details{
To access the tutorial document for this package, type in R: \code{vignette("mvGST")}
}

\value{
A list of matrices. Each matrix has possible profiles as the row names and contrasts as the column names. Ones in the appropriate cells showing which profile the gene set fit for each contrast and zeroes elsewhere. The names of the list are the names, or ID's, provided.
}
\references{
Stevens, J. R., and Isom, S. C., 2012. "Gene set testing to characterize multivariately differentially expressed genes." Conference on Applied Statistics in Agriculture
Proceedings, 24, pp. 125-137.
http://newprairiepress.org/agstatconference/2012/proceedings/10/

Mecham, D. S. (2014) "mvGST: Multivariate and Directional Gene Set Testing". MS Project, Utah State University, Department of Mathematics and Statistics.  http://digitalcommons.usu.edu/gradreports/382/
}
\author{
John R. Stevens and Dennis S. Mecham
}

\examples{
data(mvGSTsamples)

# object obatoclax.mvGST returned by profileTable

# Returns a list of 2 matrices: one matrix showing the significance profile for 
# each contrast for the gene set GO:0000003, and one matrix for GO:0000019
go2Profile(c("GO:0000003", "GO:0000019"), obatoclax.mvGST)

# See package vignette for larger examples with discussion: 
#    vignette("mvGST")
}

