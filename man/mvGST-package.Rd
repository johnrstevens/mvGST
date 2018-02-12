\name{mvGST-package}
\alias{mvGST-package}
\alias{mvGST}
\docType{package}
\title{
Multivariate and directional gene set testing
}
\description{
mvGST provides platform-independent tools to identify GO terms (gene sets) 
that are differentially active (up or down) in multiple contrasts of interest.  
Given a matrix of one-sided p-values (rows for genes, columns for contrasts), 
mvGST uses meta-analytic methods to combine p-values for all genes annotated 
to each gene set, and then classify each gene set as being significantly more 
active (1), less active (-1), or not significantly differentially active (0) 
in each contrast of interest.  With multiple contrasts of interest, each gene set 
is assigned to a profile (across contrasts) of differential activity.  
Tools are also provided for visualizing (in a GO graph) the gene sets 
classified to a given profile.
}
\details{
\tabular{ll}{
Package: \tab mvGST\cr
Type: \tab Package\cr
Version: \tab 0.1.14\cr
Date: \tab 2018-01-25\cr
License: \tab GPL-3\cr

To access the tutorial document for this package, type in R: \code{vignette("mvGST")}

}
User must provide a matrix a p-values with rows representing genes and columns representing contrasts that were tested. The contrasts must be given in the form var(1).var(2)...var(n-1).var(n) var(1) is the variable that will define the possible significance profiles. Each profile is a set of zeros, ones, and negative ones meaning significantly greater than, less than, or not significant. For example, if var(1) has levels a and b, the profile 1,0 would indicate significance (greater than) at level a and no significance at level b. The main result of \code{profileTable} is a matrix with significance profiles for row names and contrasts tested (not including var(1)) for column names and the total number of gene sets that fit each profile for each contrast in the cells.

\code{pickOut} returns the Gene Ontology ID's for the gene sets in a given cell of the \code{results.table} produced by \code{profileTable}.

\code{graphCell} uses the gene sets from \code{pickOut} to make a GO graph and displays the names of the gene sets in a legend.

\code{go2Profile} returns a matrix similar to the \code{results.table} for each desired gene set. The only difference is that there is only one gene set included in these matrices.
}
\author{
John R. Stevens, Dennis S. Mecham and Garrett Saunders

Maintainer: John R. Stevens <john.r.stevens@usu.edu>

}




\references{
Stevens, J. R., and Isom, S. C., 2012. "Gene set testing to characterize multivariately differentially expressed genes." Conference on Applied Statistics in Agriculture
Proceedings, 24, pp. 125-137.
http://newprairiepress.org/agstatconference/2012/proceedings/10/

Mecham, D. S. (2014) "mvGST: Multivariate and Directional Gene Set Testing". MS Project, Utah State University, Department of Mathematics and Statistics.  http://digitalcommons.usu.edu/gradreports/382/

Saunders, G., 2014. "Family-wise error rate control in QTL mapping and 
gene ontology graphs with remarks on family selection." PhD thesis, 
Utah State University, Department of Mathematics and Statistics.  
http://digitalcommons.usu.edu/etd/2164/
}
\keyword{ mvGST }

\examples{
# See examples in help files of four main functions:
#   ?profileTable; ?pickOut; ?go2Profile; ?graphCell

# See package vignette for larger examples with discussion: 
#    vignette("mvGST")
}
