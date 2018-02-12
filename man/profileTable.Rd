\name{profileTable}
\alias{profileTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a table of significance profiles of gene sets from given matrix of p-values
}
\description{
\code{profileTable} takes a matrix of one-sided p-values with rows representing genes and columns representing contrasts. Rows (and p-values) are combined using Stouffer's method so that the new rows represent gene sets. P-values are then adjusted for multiple hypothesis testing using the Benjamini-Yekutieli adjustment and converted to 1 (p-value < alpha/2), -1 (alpha > 1-alpha/2) or 0 (not significant). Then each gene set is classified according to its significance profile (across one of the factors) for each of the remaining contrasts.
}
\usage{
profileTable(pvals, gene.names, contrasts, list.groups, sig.level = 0.05,
             gene.ID, organism, affy.chip, ontology = "BP", method = 2, 
             minsize = 1, maxsize = Inf, mult.adj = "BY")
}

\arguments{
  \item{pvals}{
A matrix containing the one-sided p-values corresponding to the various genes (rows) and contrasts (columns)
}
  \item{gene.names}{
A character vector containing the gene names that correspond to the rows of the matrix of p-values. If NULL, then \code{gene.names <- rownames(pvals)}
}
  \item{contrasts}{
A character vector containing the contrasts that correspond to each column in the matirx of p-values. Must either be in format: Var1 or Var1.Var2 (Var2 is optional). The number of levels in Var1 determines the dimensions of the profiles (i.e. if Var1 has 2 levels, then there will be two columns for the profiles in the returned table). Var2 determines the number of columns, or strata, that will be reported in the returned table for each profile. If Var2 is not given, then there will only be one column reported, which will be the ontology chosen. If NULL, then \code{contrasts <- colnames(pvals)}
}

  \item{list.groups}{
An optional list containing user-defined gene sets
}
  \item{sig.level}{
The alpha level that should be used. Default is .05. This level is divided equally between the two sides of the test. So, for the default level of .05, any p-value less than .025 or greater than .975 is considered significant.
}
  \item{gene.ID}{
Gene naming system used for the gene names. Used to generate list of gene sets mapping genes to Gene Ontology sets. \code{gene.ID} can be "entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene", or "affy" among others.  (See \code{ID} argument in the \code{annFUN.org} function of the \code{topGO} package for supported levels.) If ID is all numeric and is not listed above, see http://biit.cs.ut.ee/gprofiler/gconvert.cgi for proper input.
}
  \item{organism}{
The organism that the genes come from. Used to generate list of gene sets mapping genes to Gene Ontology sets. The organism name should be the first letter of the scientific name and the second word of the scientific name, all lower case. For example, human would be "hsapiens".
organism must be one of the following: agambiae, athaliana, btaurus, celegans, cfamiliaris, dmelanogaster, drerio, ecoliK12, ecoliSakai, ggallus, hsapiens, mmusculus, mmulatta, pfalciparum, ptrogldytes, rnorvegicus, scerevisiae, scoelicolor, sscrofa, tgondii, or xlaevis
}
  \item{affy.chip}{
The type of affy.chip used, if gene.ID == "affy". If \code{abatch} is an AffyBatch object, then use \code{affy.chip = annotation(abatch)}
}
  \item{ontology}{
The ontology that should be used for gene sets: "BP", "MF", or "CC". Default is "BP" 
}
  \item{method}{
The method for handling gene name translation issues. Default is 2 (See "Details" below).

}
  \item{minsize}{
The minimum number of genes a gene set can have and still be included in the list of gene sets,
if the list is not provided by the user.
}
  \item{maxsize}{
The maximum number of genes a gene set can have and still be included in the list of gene sets,
if the list is not provided by the user.
}
  \item{mult.adj}{
The type of multiple hypothesis adjustment to make. BY is a Benjamini-Yekutieli adjustment. SFL
is a Short Focus Level adjustment. 
  }
}
\details{
User must provide a matrix a p-values with rows representing genes and columns representing contrasts that were tested. The contrasts must be given in the form Var1.Var2 or just Var1 (Var2 is optional). The possible significance profiles will be defined by Var1. The number of dimensions of the profiles is the same as the number of levels of Var1. Each profile is a set of zeros, ones, and negative ones meaning significantly greater than, less than, or not significant. For example, if Var1 has levels a and b, the profile 1,0 would indicate significance (greater than) at level a and no significance at level b. The main result of \code{profileTable} is a matrix with significance profiles for row names and contrasts tested (not including Var1) for column names and the total number of gene sets that fit each profile for each contrast in the cells.

If the gene names given are not affy or entrez ID's, then they will have to be translated to entrez ID's if the user does not provide their own list of gene sets. Translation is done using \code{gconvert} from the \code{gProfileR} package. This may lead to one gene being translated to many, or many being translated to one. These problems are handled using a method of the user's choice (See Section 2.1.1 of Mecham, 2014):  
  Method 1 does nothing. As a result, some rows of p-values will be duplicated when one name translates to many. Some rows will also have the same gene name when many names translate to just one.
  Method 2 uses Stouffer's inverse normal method to combine p-values when many names translate to just one. 
  Method 3 accounts for when one name translates to many. Instead of duplicating rows of p-values, only the first of the new names is used.
  Method 4 combines methods 2 and 3. First method 2 is performed, then method 3. 
 
To access the tutorial document for this package, type in R: \code{vignette("mvGST")}
}

\section{warning}{
If user does need to use gene ID's other than affy or entrez, it is strongly recommended that the user provides his or her own list of gene sets. Some of the translations that have been tested can take a VERY long time.
}

\value{
Returns an object of class "mvGST". An object of class "mvGST" is a list containing the following components:
\item{results.table}{
A matrix with possible profiles as row names and contrasts as column names. The cells of the matrix show how many gene sets have each profile for each contrast.
}
\item{raw.pvals}{
A matrix of the original p-values provided with gene names as row names and contrasts as column names.
}
\item{grouped.raw}{
A matrix of Stouffer combined p-values. Each row is for a gene set and each column is for a contrast. 
}
\item{adjusted.group.pvals}{
The same matrix as in grouped.raw, but with a the requested adjustment ("BY" or "SFL") being performed within each column.
}
\item{ones.zeroes}{
A matrix showing the significance results of each of the BY adjusted p-values. 1 means significantly greater. -1 mean significantly less. 0 means not significant.
}

\item{ord.lev}{
The levels of the ordered  variable (the variable that defines the profiles).
}
\item{contrasts}{
The contrasts that are the column names of $results.table.
}
\item{group.names}{
The Gene Ontology ID's of the gene sets in the selected ontology.
}

}
\references{
Stevens, J. R., and Isom, S. C., 2012. "Gene set testing to characterize multivariately differentially expressed genes." Conference on Applied Statistics in Agriculture
Proceedings, 24, pp. 125-137.
http://newprairiepress.org/agstatconference/2012/proceedings/10/

Mecham, D. S. (2014) "mvGST: Multivariate and Directional Gene Set Testing". MS Project, Utah State University, Department of Mathematics and Statistics.  http://digitalcommons.usu.edu/gradreports/382/

Saunders, G., 2014. "Family-wise error rate control in QTL mapping and gene ontology graphs with remarks on family selection." PhD thesis, Utah State University, Department of Mathematics and Statistics.  http://digitalcommons.usu.edu/etd/2164/
}
\author{
John R. Stevens and Dennis S. Mecham
}

\examples{
# A matrix of p-values with 3 rows (one for each gene) and 4 columns (one for
# each contrast)
pvals <- matrix(c(.001, .03, .7, .5, .01, .0002, .01, .85, .97, .99, 1, .98), 
                nrow = 3)

# A character vector with the IDs of the genes (in this case, 3 Entrez IDs)
gene.names <- c("8510", "4361", "10111")
rownames(pvals) <- gene.names

# A character vector with the contrasts that were tested (the first factor
# had levels 1 and 2; the second factor had levels a and b; in each of the
# first factor's levels, second factor levels a and b were tested 
# against a third level)
contrasts <- c("1.a", "1.b", "2.a", "2.b")
colnames(pvals) <- contrasts

# Creates an mvGST object with a results.table showing how many 
# biological process (BP) gene sets fit each significance profile 
# for each contrast.  Gene names are from Entrez, and data are human.
test <- profileTable(pvals, gene.ID = "entrez", organism = "hsapiens", 
                      ontology = "BP")

# See how many gene sets were classified into each profile
# (across levels of the first factor) for each level of the
# second factor.  For example, the '0  0' line represents the 
# profile of gene sets that are not differentially active in either
# level of the first factor.
test
					  
# See package vignette for larger examples with discussion: 
#    vignette("mvGST")
}

