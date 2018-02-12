\name{pickOut}
\alias{pickOut}

\title{
Returns gene sets in a cell of the \code{results.table}.
}
\description{
\code{pickOut} returns a character vector with the Gene Ontology ID's of the gene sets with a particular significance profile for a certain contrast (the gene sets in one cell of the \code{results.table} of an mvGST object).
}
\usage{
pickOut(mvgst, row, col = 1)
}

\arguments{
  \item{mvgst}{
A mvGST object. \code{mvGST$results.table} must not be NULL.
}
  \item{row}{
The row of the desired profile.
}
  \item{col}{
The column of the desire contrast. Column refers to the levels of Var2, if Var2 was used. It is the number of the column after the 1, 0, -1 columns that show the profiles. Default value is 1. 
}
}

\details{
To access the tutorial document for this package, type in R: \code{vignette("mvGST")}
}

\value{
A data frame containing the ID's of the gene sets in the given row and column along with the GO descriptions of each gene set and the adjusted p-values for each contrast tested.
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
obatoclax.mvGST

# returns the GO ID's from the cell in the third row and first column\
# (the column for cell line RS4) of the results.table of 
# the object returned by profileTable
pickOut(obatoclax.mvGST, 3, 1)

# See package vignette for larger examples with discussion: 
#    vignette("mvGST")

}

