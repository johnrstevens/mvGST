\name{graphCell}
\alias{graphCell}

\title{
Makes a GO graph highlighting the GO terms in the selected cell of the results table.
}
\description{
\code{graphCell} relies on tools in the Rgraphviz package.
\code{graphCell} uses \code{pickOut} to get the GO terms in a specific cell of the results table. A GO graph is created from those GO terms, and can be interactive if desired. Also, if desired, a legend showing the names of the GO terms can be printed. If the graph is interactive, use esc to end interaction with graph.
}
\usage{
graphCell(object, row, col = 1, set, ontology = "BP", interact = TRUE, 
          legend.pos = "bottomleft", print.legend = TRUE,
		  use.col="red", bg.col = "grey80")
}

\arguments{
  \item{object}{
A \code{mvGST} object with a final \code{results.table}
}
  \item{row}{
The row of the desired cell.
}
  \item{col}{
The column of the desired cell. Column refers to the levels of Var2, if Var2 was used. It is the number of the column after the 1, 0, -1 columns that show the profiles. Default value is 1. 
}
  \item{set}{
Optional argument that is a data frame with the first column containing the GO ID's that should be used to make the GO graph. The data frame returned by \code{pickOut} can be used.  
}
  \item{ontology}{
The ontology, within Gene Ontology, that should be used ("BP", "MF", "CC").
}
  \item{interact}{
Indicates whether or not the graph should be interactive. If interactive, use esc to end interaction with graph.
}
  \item{legend.pos}{
If interactive, indicates the desired position of the legend that shows name and GO ID of selected node.
}
  \item{print.legend}{
Indicates if the legend should also be printed separately, showing GO names of all nodes.
}
  \item{use.col}{
Color to highlight the nodes representing gene sets of interest in the resulting graph.
}
  \item{bg.col}{
Color to use for the "background" in the graph when focusing on the gene sets of interest.  This is the color used for the border of all nodes, the labels of all nodes NOT representing gene sets of interest, and all edges.
}
}

\value{
  Invisibly returns \code{NULL}.
}

\details{
To access the tutorial document for this package, type in R: \code{vignette("mvGST")}
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

# plots a GO Graph highlighting the GO ID's from the cell 
# in the fifth row and first column (the column for 
# cell line RS4) of the results.table of the object 
# returned by profileTable
graphCell(obatoclax.mvGST, 5, 1, ontology = "BP", interact = FALSE)

# See package vignette for larger examples with discussion: 
#    vignette("mvGST")
}

