\name{mvGST.other}
\alias{changeTO10}
\alias{combinePvalues}
\alias{convertPvalues}
\alias{cut}
\alias{fillInList}
\alias{finalResults}
\alias{geneNameConvertRows}
\alias{generateGeneSets}
\alias{go2GeneSet}
\alias{hartung}
\alias{interactiveGraph}
\alias{method2}
\alias{method3}
\alias{method4}
\alias{mvGSTObject}
\alias{mvSort}
\alias{oneSideBYAdjust}
\alias{print.mvGST}
\alias{print.summary.mvGST}
\alias{profileCombine}
\alias{profiles}
\alias{separate}
\alias{summary.mvGST}
\alias{tableColumns}
\alias{distributeWeight}
\alias{getCurrentChildren}
\alias{makeCoherent}
\alias{getAncestorsAndOffspring}
\alias{turnListAround}

\title{
mvGST Other Functions
}
\description{
Internal functions used by the main functions of \code{mvGST} (\code{graphCell}, \code{profileTable}, \code{go2GeneSets}, \code{pickOut}, and \code{p.adjust.SFL}):
\tabular{ll}{
\tab \cr
\code{changeTO10}  \tab changes p-values in matrix to -1, 0, or 1, indicating significance. \cr
\code{combinePvalues}  \tab Uses Stouffer's method to combine p-values in gene sets \cr
\code{convertPvalues}  \tab converts p-values from one-sided to two-sided and vice versa \cr
\code{cut}  \tab cuts out any all 0 rows from the matrix of results \cr
\code{fillInList}  \tab fills in list of gene sets, ensuring that any gene in a child set is also in the parent set \cr
\code{finalResults}  \tab counts the number of gene sets corresponding to each profile and contrast and creates \code{results.table} \cr
\code{geneNameConvertRows}  \tab handles one-to-many and many-to-one problems of gene name translation \cr
\code{generateGeneSets}  \tab creates list of gene sets \cr 
\code{go2GeneSet}  \tab creates a single matrix showing which profile a gene set fits into for each contrast tested \cr
\code{hartung}  \tab Combines a set of p-values using Hartung's modified inverse normal method \cr
\code{interactiveGraph}  \tab Creates a GO Graph of selected GO terms \cr
\code{method2}  \tab Uses Hartung's method to combine p-values and handles many-to-one problem of gene name translation \cr
\code{method3}   \tab Handles one-to-many problem of gene name translation by eliminating all but one of the many\cr
\code{method4}   \tab Combines \code{method2} and \code{method3} \cr
\code{mvGSTObject}   \tab Creates an mvGST object \cr
\code{mvSort}  \tab Sorts \code{results.table} by row totals \cr
\code{oneSideBYAdjust}  \tab Convert p-values to two-sided, then uses Benjamini-Yekutieli adjustment, and converts back to one-sided \cr
\code{print.mvGST}    \tab prints the \code{results.table} in a clean way \cr
\code{print.summary.mvGST}    \tab prints the summary of an mvGST object \cr
\code{profileCombine}   \tab Determines profiles of each gene set at each contrast \cr
\code{profiles}  \tab Defines the possible profiles \cr
\code{separate}   \tab splits genes (and their p-values) into gene sets \cr
\code{summary.mvGST}    \tab creates a summary of an mvGST object \cr
\code{tableColumns}   \tab Creates the column names for the \code{results.table} \cr
\tab \cr
\code{distributeWeight} \tab helper function for p.adjust.SFL \cr
\code{getCurrentChildren} \tab helper function for p.adjust.SFL \cr
\code{makeCoherent} \tab helper function for p.adjust.SFL \cr
\code{getAncestorsAndOffspring} \tab helper function for p.adjust.SFL \cr
\code{turnListAround} \tab From globaltest package, helper function for p.adjust.SFL \cr
  }
}

\details{
To access the tutorial document for this package, type in R: \code{vignette("mvGST")}
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
\author{
John R. Stevens, Dennis S. Mecham, Garrett Saunders
}






