separate <- function(pmat, list.groups){
  # Uses Stouffer's method to combine p-values for gene sets
  #
  # Args:
  #   pmat: An MVGST object with a matrix of p-values with corresponding 
  #         gene names as the row names
  #   list.groups: A list of character vectors containing the gene sets
  #
  # Returns:
  #   An MVGST object with a matrix of p-values with corresponding gene 
  #   set names st the row names
  mvgst <- pmat$raw.pvals
  new.ps <- matrix(NA, nrow = length(list.groups), ncol = ncol(mvgst))
  for (i in seq_along(list.groups)){
    t <- is.element(rownames(mvgst), list.groups[[i]])
    if (any(t)){
      new.ps[i, ]<- apply(matrix(mvgst[t, ], nrow = sum(t)), MARGIN = 2, FUN = combinePvalues)
    } 
  }

  new.ps <- matrix(new.ps, dimnames = list(rep(1:nrow(new.ps)),
                                           colnames(mvgst)),
                   nrow = nrow(new.ps))
  
  pmat$grouped.raw <- new.ps
  return(pmat)
}


####################################################

finalResults <- function(pmat){
  # Takes a matrix of significance results (-1, 0, 1) and creates a table
  # with the total number of gene sets that fall into each possible 
  # profile
  #
  # Args:
  #   pmat: An MVGST object containing the matrix of ones, zeroes, and
  #         negative ones indicating whether each gene set is
  #         significantly less differentially expressed, not 
  #         significantly differentially expressed, or significantly
  #         greater differentially expressed for each contrast
  #
  # Returns:
  #   An MVGST object containing a matrix with row names that are each
  #   possible significance profile for the first variable, column names 
  #   that are each combination of the remaining variables, and cell 
  #   values that indicate the number of gene sets that have the 
  #   corresponding significance profile. 
  mat <- pmat$ones.zeroes
  # setting up the final results matrix
  profile <- profiles(mat)
  columns <- tableColumns(mat)
  observed <- profileCombine(mat)
  
  profs <- as.list(rep(NA, nrow(profile)))
  # counting the number of gene sets that fit in each cell of matrix
  nrowprof <- nrow(profile)
  for (i in seq_len(nrowprof)){
    profs[[i]] <- profile[i,]
  }
  final <- matrix(NA, nrow = nrow(profile), ncol = length(columns),
                  dimnames = list(profs, columns))
   nrfinal <- nrow(final)
   ncfinal <- ncol(final)   
  for (i in seq_len(nrfinal)){
    for (j in seq_len(ncfinal)){
      final[i,j] <- sum(apply(observed[[j]], MARGIN = 1, FUN = function(x) all(x == profile[i,])), na.rm = TRUE)
    }
  }
  pmat$results.table <- final
  pmat$ord.lev <- dimnames(profile)[[2]]
  pmat$contrasts <- dimnames(final)[[2]]
  return(pmat)
}

###############################################################

profiles<-function(mat){
  # Takes a matrix of significance results (-1, 0, 1) and creates a matrix
  # with each row being a possible significance profile and the number of
  # rows being equal to the number or possible significance profiles.
  #
  # Args:
  #   mat: A matrix with column names that are the in the format 
  #        Var1.Var2.Var3...Var(n-1).Var(n) that represent the contrasts
  #        being tested.
  #
  # Returns:
  #   A matrix with each row being a possible significance profile and the 
  #   number of rows being equal to the number or possible significance 
  #   profiles. 
  names <- dimnames(mat)[[2]]
  prof <- rep(NA, length(names))
  for (i in seq_along(names)){
    prof[i] <- unlist(strsplit(names[i], "[.]"))[1]
  }
  t <- duplicated(prof)
  unique <- prof[!t]
  nprofiles <- 3 ^ length(unique)
  profiles <- matrix(rep(NA), ncol = length(unique),
                     nrow = nprofiles, 
                     dimnames = list(rep(1:nprofiles), unique))
  for (i in seq_along(unique)){
    profiles[, i] <- rep(c(rep(-1, nprofiles / 3 ^ i),
                           rep(0, nprofiles / 3 ^ i),
                           rep(1, nprofiles / 3 ^ i)), 3 ^ (i - 1))
  }
  return(profiles)
}

#################################################################


mvGSTObject <- function(gene.names, contrasts, pvals, groups){
  # Creates an object of class mvGST
  #
  # Args:
  #   gene.names: A character vector containing the gene names that correspond
  #               to the rows of the matrix of p-values
  #   contrasts: A character vector containing the contrasts that correspond to
  #              each column in the matirx of p-values. Must be in format:
  #              Var1.Var2.Var3...Var(n-1).Var(n)
  #   pvals: A matrix containing the p-values corresponding to the various genes
  #          and contrasts
  #   groups: An optional list containing user-defined gene sets 
  #
  # Returns:
  #   A mvGST object with most things still empty
  object <- matrix(pvals, dimnames = list(gene.names, contrasts),
                   nrow = length(gene.names))
  obj <- list(raw.pvals = object,
              results.table = matrix(rep(NA,4),nrow=2),
              ord.lev = NA,
              contrasts = NA,
              grouped.raw = NA,
              adjusted.group.pvals = NA,
              ones.zeroes = NA,
              group.names = names(groups))
  class(obj) <- "mvGST"
  return(obj)
}

######################################################################

tableColumns <- function(mat){
  # Takes a matrix of significance results (-1, 0, 1) and returns the 
  # contrasts that need to be the column names for the final results
  # matrix
  #
  # Args:
  #   mat: A matrix with column names that are the in the format 
  #        Var1.Var2 or Var1that represent the contrasts
  #        being tested.
  #
  # Returns:
  #   A character vector in the format Var2 or "ontology"
  names <- dimnames(mat)[[2]]
  col <- rep(NA, length(names))
  columns <- rep(NA, length(names))
  for (i in seq_along(names)){
    col[i] <- unlist(strsplit(names[i], "[.]"))[1]
    columns[i] <- substr(names[i], nchar(col[i]) + 2, nchar(names[i]))
  }
  t <- duplicated(columns)
  unique <- columns[!t]
  return(unique)
}

##########################################################################

convertPvalues <- function(pvals, relativity, two.sided = TRUE){
  # Converts a matrix of p-values from two-sided to one-sided or vice versa.
  #
  # Args:
  #   pvals: A matrix of p-values.
  #   relativity: Only used when two.sided == TRUE. Numeric value that is 
  #               greater than 0 if the one-sided p-value should be less 
  #               than .5.
  #   two.sided: TRUE if pvals contains two-sided pvalues, FALSE if pvals
  #              contains one-sided pvalues.

  if(two.sided){
    newP <- ifelse(relativity < 0, 1 - pvals / 2, pvals / 2)
  } else { 
    newP <- ifelse(pvals < .5, pvals * 2, (1 - pvals) * 2)
  }
  return(newP)
}

################################################################

combinePvalues <- function(pvals){
  # Uses Stouffer's method to combine p-values
  #
  # Args:
  #   pvals: A vector of p-values.
  #
  # Returns:
  #   A single combined p-value
  comb.p <- pnorm(sum(qnorm(pvals)) / sqrt(length(pvals)))
  return(comb.p)
}

##################################################################

profileCombine <- function(mat){
  # Takes a matrix of significance results (-1, 0, 1) and list of matrices 
  # with each matrix containing the significance profiles for a given
  # Var2.Var3...Var(n-1).Var(n)
  #
  # Args:
  #   mat: A matrix with column names that are the in the format 
  #        Var1.Var2.Var3...Var(n-1).Var(n) that represent the contrasts
  #        being tested.
  #
  # Returns:
  #   list of matrices with each matrix containing the significance 
  #   profiles for a given Var2.Var3...Var(n-1).Var(n)
  names <- dimnames(mat)[[2]]
  col <- rep(NA, length(names))
  duplicates <- rep(NA, length(names))
  columns <- tableColumns(mat)
  each <- profiles(mat)
  mats <- rep(list(matrix(rep(NA), nrow = nrow(mat),
                          ncol = ncol(each))), length(columns))
  for(i in seq_along(columns)){
    for (j in seq_along(names)){
      col[j] <- unlist(strsplit(names[j], "[.]"))[1]
      duplicates[j] <- substr(names[j], nchar(col[j]) + 2, nchar(names[j]))
    }
    t <- duplicates == columns[i]
    mats[[i]] <- matrix(mat[, t], nrow = nrow(mat))
  }
  return(mats)
}

############################################################

changeTO10 <- function(pvals.mat, sig.level){
  # Converts a matrix of one-sided p-values to matrix of ones, zeroes and
  # negative ones where 1 represents significantly greater, 0 represents
  # no significance, and -1 represents significantly less than
  #
  # Args:
  #   pvals.mat: A mvGST object containing a matrix of one-sided p-values.
  #   sig.level: Significance level to be used. P-values less than 
  #              sig.level / 2 are converted to 1. P-values greater
  #              than 1 - sig.level / 2 are converted to -1.
  #
  # Returns:
  #   A matrix of ones, zeroes, and negative ones.
  pvals <- pvals.mat$adjusted.group.pvals
  ind.profile <- matrix(rep(0), nrow = nrow(pvals), ncol = ncol(pvals),
                        dimnames = dimnames(pvals))
  ind.profile <- ifelse(pvals <= sig.level / 2, 1, 
                        ifelse(pvals >= 1 - sig.level / 2, -1, 0))
  pvals.mat$ones.zeroes <- ind.profile
  return(pvals.mat)
}

#########################################################################

oneSideBYAdjust<-function(pval.mat){
  # Converts one-sided p-values to two-sided. Performs a Benjamini-Yekutieli
  # adjustment on the two-sided p-values. Converts the adjusted two-sided 
  # p-values back to one-sided p-values.
  #
  # Args:
  #   pval.mat: A mvGST object that contains a matrix of Stouffer combined 
  #             p-values
  #
  # Returns:
  #   A mvGST object that also contains a matrix of BY adjusted, Stouffer 
  #   combined p-values.
  pvals <- pval.mat$grouped.raw

  two.sided <- convertPvalues(pvals, two.sided = FALSE)
  two.adjusted <- apply(two.sided, MARGIN = 2, FUN = p.adjust, method = "BY")
  relative <- ifelse(pvals < .5, 1, -1)
  one.combined <- convertPvalues(two.adjusted, relative)
  pval.mat$adjusted.group.pvals <- matrix(one.combined, nrow = nrow(pvals),
                                          dimnames = dimnames(pvals))
  return(pval.mat)
}

################################################################

cut <- function(y){
  # Removes the rows with all zeroes from the final results table in a mvGST object.
  #
  # Args:
  #   y: A mvGST object that contains a final results matrix
  #
  # Returns:
  #   A mvGST object that contains a final results matrix with no all zero rows
  x <- y$results.table
  temp.mat <- matrix(0, ncol = ncol(x))
  temp.names <- NA
  nrx <- nrow(x)
  for (i in seq_len(nrx)){
    if (sum(x[i,]) != 0){
      temp.mat <- rbind(temp.mat, x[i, ])
      temp.names <- c(temp.names, dimnames(x)[[1]][i])
    }
  }
  temp.mat <- temp.mat[-1, ]
  temp.names <- temp.names[-1]
  final.mat <- matrix(temp.mat, ncol = ncol(x), 
                      dimnames = c(list(temp.names), list(dimnames(x)[[2]])))
  y$results.table <- final.mat
  return(y)
}


#####################################################################

mvSort <- function(y){
  # Sorts the rows in the final results table in a mvGST object from the greatest
  # row total to the least.
  #
  # Args:
  #   y: A mvGST object that contains a final results matrix
  #
  # Returns:
  #   A mvGST object that contains a final results matrix with rows sorted by
  #   row total.
  x <- y$results.table
  temp <- x
  row.sum <- apply(temp, MARGIN = 1, FUN = sum)
  ranked <- rank(1 / row.sum, ties.method = "first")
  nrt <- nrow(temp)
  for (i in seq_len(nrt)){
    x[ranked[i], ] <- temp[i, ]
    dimnames(x)[[1]][ranked[i]] <- dimnames(temp)[[1]][i]
  }
  y$results.table <- x
  return(y)
}

#########################################################################

print.mvGST <- function(x, ...){
  # Prints an object of class mvGST
  #
  # Args:
  #   x: A mvGST object 
  
  #!# Allow for all-NA profile (like if a GO id passed to go2Profile was not among those actually tested)
  if(dim(x$results.table)[1]==0){
      NAprof <- matrix(NA, nrow=1, ncol=length(x$ord.lev))
	  colnames(NAprof) <- x$ord.lev
	  fillcols <- matrix(1,nrow=1, ncol=dim(x$results.table)[2])
      colnames(fillcols) <- colnames(x$results.table)
	  outNAprof <- cbind(NAprof, fillcols)
      print(outNAprof)
      warning('Gene set name(s) requested not among those tested')	  

  } else if (substr(rownames(x$results.table)[1], 1, 1) != "c"){
    print(x$results.table)
   } else {
    y <- x$results.table
    col.names1 <- dimnames(y)[[2]]
    col.names <- ifelse(nchar(col.names1) < 4, paste(col.names1, "  "), col.names1)
    row.spaces <- nchar(col.names)
    row.names <- dimnames(y)[[1]]
    col.spaces <- max(nchar(row.names))
    for (i in seq_along(row.names)){
      temp.name <- substr(row.names[i], 3, nchar(row.names[i])-1)
      row.names[i] <- " "	
      counter <- 1
	  nct <- nchar(temp.name)
      for (j in seq_len(nct)){
       char <- substr(temp.name, j, j)
       if (char == "-"){
         row.names[i] <- paste(row.names[i], char, sep = "")
       } else {
         if (char == "0"){
           row.names[i] <- paste(row.names[i], char)
		   ncxo <- nchar(x$ord.lev[counter])
           for (k in seq_len(ncxo)){
             row.names[i] <- paste(row.names[i], "")
             
           }
           counter <- counter + 1
         } else {
           if (char == "1"){
             if (substr(temp.name, j - 1, j - 1) == "-"){
               row.names[i] <- paste(row.names[i], char, sep = "")
             } else {
               row.names[i] <- paste(row.names[i], char)
             }
			 ncxol <- nchar(x$ord.lev[counter])
             for (k in seq_len(ncxol)){
               row.names[i] <- paste(row.names[i], "")
               
             }
             counter <- counter + 1
           }
         }
       }
      }
    }
    row.names <- gsub(",", " ", row.names)
    col.spaces <- max(nchar(row.names))
    cat(" ")
    for (i in seq_along(x$ord.lev)){
      cat(x$ord.lev[i], " ")
    }

    cat(col.names, "\n")
    for (i in seq_along(row.names)){
      cat(row.names[i], "")
	  ncy <- ncol(y)
      for (j in seq_len(ncy)){
        spaces <- 1 + row.spaces[j] - nchar(y[i, j]) 
        cat(y[i, j], rep("", max(spaces, 1)))
      }
      cat("\n")
    }
  }
}

##################################################################

summary.mvGST <-
  function(object, ...){
    # Creates a summary of the mvGST object, giving number of gene sets tested,
    # levels of the ordered factor, number of other factors, number of possible
    # profiles, number of profiles that have gene sets, and number of contrasts
    # tested.
    #
    # Args:
    #   object: A mvGST object
    y <- object$results.table
    gene.sets <- sum(y[, 1])
    ordered.vars <- length(object$ord.lev)
    summ <- c(gene.sets, 3 ^ ordered.vars, nrow(y), ncol(y))
    class(summ) <- "summary.mvGST"
    return(summ)
  }
###################################################

print.summary.mvGST <- function(x, ...){
  # Prints an object of class summary.mvGST 
  #
  # Args:
  #   x: A summary.mvGST object
  cat("", x[1], "gene sets \n",  x[4], "possible profiles \n", 
      x[5], "profiles used \n", x[6], "strata \n")
}

#######################################################

geneNameConvertRows <- function(pvals, gene.names, new.names, method = 1){
  # Handles translation from one set of gene names to another. 
  #
  # Args:
  #   pvals: A matrix of p-values.
  #   gene.names: A character vecter containing the original gene names.
  #   new.names: A matrix with the first column containing the original gene
  #              names, and the second column containing the corresponding new
  #              names. 
  #   method: A number from 1 to 4 that indicates what method should be used to 
  #           handle duplicates in the name translation.
  #           Method 1 does nothing. As a result, some rows of p-values will be 
  #             duplicated when one name translates to many. Some rows will also
  #             have the same gene name when many names translate to just one.
  #           Method 2 uses Hartung's modified inverse normal method to combine
  #             p-values when many names translate to just one.
  #           Method 3 accounts for when one name translates to many. Instead of 
  #             duplicating rows of p-values, only the first of the new names is
  #             used.
  #           Method 4 combines methods 2 and 3. First method 2 is performed, then
  #             method 3.
  #
  # Returns:
  #   A list containing the new p-value matrix, the new gene names, and the old
  #   gene names if method != 1.
  
  if (method == 1){
    new.pvals <- matrix(NA, nrow = nrow(new.names), ncol = ncol(pvals))
	nrnn <- nrow(new.names)
    for (i in seq_len(nrnn)){
      new.pvals[i, ] <- pvals[new.names[i, 1] == gene.names]
    }
    return(list(p.mat = new.pvals, genes = new.names[, 2]))
  }
  if (method == 3){
    return(method3(pvals, gene.names, new.names))  
  }
  if (method == 2){
    return(method2(pvals, gene.names, new.names))
  }
  if (method == 4){    
    method.two <- method2(pvals, gene.names, new.names)
    new.new.names <- data.frame(method.two$old.names, method.two$genes)
    return(method4(method.two$p.mat, method.two$old.names, new.new.names))
  }
}

###########################################################

hartung <- function(pvals){
  # Uses Hartung's modified inverse normal method to combine a set of 
  # p-values
  #
  # Args:
  #   pvals: A vector of p-values
  #
  # Returns:
  #   A single p-value.
  n <- length(pvals)
  t <- qnorm(pvals)
  rho.hat <- 1 - sum((t - sum(t) / n) ^ 2) / (n - 1)
  rho.hat.star <- max(-1 / (n - 1), rho.hat)
  kappa <- .1 * (1 + 1 / (n-1) - rho.hat.star)
  combined.test.statistic <- sum(t) / 
    sqrt(n + (n ^ 2 - n) * 
           (rho.hat.star + kappa * sqrt(2 / (n+1)) * (1 - rho.hat.star)))
  return(pnorm(combined.test.statistic))
}

############################################################

method3 <- function(pvals, gene.names, new.names){
  # Accounts for the one-to-many gene name translation issue
  # by only using the first of the many new names that are 
  # possible.
  #
  # Args:
  #   pvals: A matrix of p-values.
  #   gene.names: A character vector with the old gene names
  #   new.names: A data frame with old gene names in the first
  #              column and corresponding new gene names in the
  #              second column.
  #
  # Returns:
  #   A list with the new matrix of p-values and the corresponding
  #   gene names.
  new.names <- new.names[order(new.names[, 1]), ]
  index <- rep(TRUE, nrow(new.names))
  for (i in seq_along(new.names[, 1])[-1]){
    if (any(new.names[i, 1] == new.names[1:(i - 1), 1])){
      index[i] <- FALSE
    }
  }
  trunc.new.names <- new.names[index, ]
  new.pvals <- matrix(NA, nrow = length(index), ncol = ncol(pvals))
  for (i in seq_along(index)){
    if (index[i]){
      new.pvals[i, ] <- pvals[new.names[i, 1] == gene.names]
    }
  }
  new.pvals <- new.pvals[!is.na(new.pvals[, 1]),]
  return(list(p.mat = new.pvals, genes = trunc.new.names[, 2]))
}

######################################################################

method2 <- function(pvals, gene.names, new.names){
  # Accounts for the many-to-one gene name translation issue
  # by using Hartung's modified inverse normal method to combine
  # the p-values
  #
  # Args:
  #   pvals: A matrix of p-values.
  #   gene.names: A character vector with the old gene names
  #   new.names: A data fram with old gene names in the first
  #              column and corresponding new gene names in the
  #              second column.
  #
  # Returns:
  #   A list with the new matrix of p-values and the corresponding
  #   gene names.
  new.names <- new.names[order(new.names[, 2]), ]
  trunc.new.names <- unique(new.names[, 2])
  matches <- rep(0, nrow(new.names))
  temp.names <- c("NOT A GENE", as.character(new.names[, 2]), "NOT A GENE")
  for(i in seq_along(new.names[, 2])[-1]){
    if (temp.names[i] == temp.names[i + 1] & temp.names[i] != temp.names[i - 1]){
      j <- 1
      k <- i
      while (temp.names[k] == temp.names[k + 1]){
        j <- j + 1
        k <- k + 1
      }
      matches[i - 1] <- j
    }
  }
  old.names <- character(length(trunc.new.names))
  new.p.mat <- matrix(NA, ncol=ncol(pvals), nrow=length(trunc.new.names))
  j <- 1
  k <- 0
  for (i in seq_along(matches)){
    if (j > i) {k <- k + 1; next}
    if (matches[i] == 0){
      new.p.mat[i-k,] <- pvals[gene.names == new.names[j, 1]]
      old.names[i-k] <- new.names[j, 1]
      j <- j + 1
    } else {
      originals <- new.names[new.names[j,2] == new.names[,2], 1]
      temp <- rep(FALSE, length(gene.names))
      for (l in seq_along(originals)){
        if (any(originals[l] == gene.names)){
          temp[originals[l] == gene.names] <- TRUE
        }
      }
      if (!is.matrix(pvals[temp, ])){
        new.p.mat[i-k,] <- combinePvalues(pvals[temp,])
      } else {
       new.p.mat[i-k,] <- apply(pvals[temp, ], MARGIN = 2, 
                               FUN = combinePvalues)
      }
      old.names[i] <- paste("other", as.character(i))
      j <- j + matches[i]
    }
  }
  old.names <- old.names[!is.na(old.names) & old.names != ""]
  return(list(p.mat = new.p.mat, genes = trunc.new.names, old.names = old.names))
}

#####################################################################################

method4 <- function(pvals, gene.names, new.names){
  # Does what method3 does, but for a p-value matrix that has already
  # gone through method2. 
  #
  # Args:
  #   pvals: A matrix of p-values.
  #   gene.names: A character vector with the old gene names
  #   new.names: A data fram with old gene names in the first
  #              column and corresponding new gene names in the
  #              second column.
  #
  # Returns:
  #   A list with the new matrix of p-values and the corresponding
  #   gene names.
  temp <- new.names
  new.names <- new.names[order(temp[,1]),]
  gene.names <- gene.names[order(temp[, 1])]
  pvals <- pvals[order(temp[, 1]), ]
  if (!is.matrix(pvals)){
    pvals <- matrix(pvals)
  }
  old.names <- c("NOT A GENE", gene.names, "NOT A GENE")
  keepers <- logical()
  for(i in seq_along(old.names[-1])[-1]){
    if (old.names[i] != old.names[i - 1] & old.names[i] != old.names[i + 1]){
      keepers[i - 1] <- TRUE
    } else {
      keepers[i - 1] <- FALSE
    }
  }
  keepers[length(keepers) + 1] <- TRUE
  k <- 0
  j <- 1
  for (i in seq_along(keepers[-1])){
    if (keepers[i] == TRUE){
      next
    } else {
      if (k >= i){
        next
      } else {
        j <- 1
        while (!keepers[i + j]){
          j <- j + 1
        }
        keepers[i] <- TRUE
        k <- i + j - 1
      }
    }
  }
  keepers <- keepers[-length(keepers)]     
  return(list(p.mat = pvals[keepers, ], genes = new.names[keepers, 2]))    
}

###########################################################################

generateGeneSets <- function(ontology, species, ID, affy.chip){
  # Generates a list of gene sets based on gene ontology.
  #
  # Args:
  #   ontology: The specific ontology within to be used. Either
  #             "BP", "MF", or "CC".
  #   species: The organism being studied. It is made up of the
  #            first letter of the scientific name and the last
  #            word of the scientific name. For example, human is
  #            "hsapien"
  #   ID: The naming system being used on the genes
  #   affy.chip: If ID = "affy", this is the specific chip that was
  #              used
  #
  # Returns:
  #   A list of character vectors. Each vector contains the names
  #   of the genes in a set. The names of the elements of the list
  #   are the Gene Ontology ID's for each gene set.
  if (ID == "affy"){
    gene.set.list <- annFUN.db(ontology, affyLib = affy.chip)
    return(gene.set.list)
  } else {
    
    ### define species info with species.db <-    ###
	species.db <- switch(species,
       agambiae      = "org.Ag.eg.db",
       athaliana     = "org.At.tair.db",
       btaurus       = "org.Bt.eg.db",
       celegans      = "org.Ce.eg.db",
       cfamiliaris   = "org.Cf.eg.db",
       dmelanogaster = "org.Dm.eg.db",
       drerio        = "org.Dr.eg.db",
       ecoliK12      = "org.EcK12.eg.db",
       ecoliSakai    = "org.EcSakai.eg.db",
       ggallus       = "org.Gg.eg.db",
       hsapiens      = "org.Hs.eg.db",
       mmusculus     = "org.Mm.eg.db",
       mmulatta      = "org.Mmu.eg.db",
       pfalciparum   = "org.Pf.plasmo.db",
       ptroglodytes  = "org.Pt.eg.db",
       rnorvegicus   = "org.Rn.eg.db",
       scerevisiae   = "org.Sc.sgd.db",
       scoelicolor   = "org.Sco.eg.db",
       sscrofa       = "org.Ss.eg.db",
       tgondii       = "org.Tgondii.eg.db",
       xlaevis       = "org.Xl.eg.db",
	   stop("species must be one of the following: 
        agambiae, athaliana, btaurus, celegans, cfamiliaris,
        dmelanogaster, drerio, ecoliK12, ecoliSakai, ggallus,
        hsapiens, mmusculus, mmulatta, pfalciparum, ptrogldytes,
        rnorvegicus, scerevisiae, scoelicolor, sscrofa,
        tgondii, or xlaevis")
	   )
    
    if (any(ID == c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene"))){
      gene.set.list <- annFUN.org(ontology, mapping = species.db, ID = ID)
      return(gene.set.list)
    } 
  }
}

################################################################

interactiveGraph <- function(GO.Graph, color.nodes, interact = FALSE,
                              legend.pos = "bottomleft", print.legend = FALSE,
							  use.col="red", bg.col = "grey80"){
  # Creates a graph showing all given gene sets and all parent gene sets
  #
  # Args: 
  #   GO.Graph: A graphNEL object created by the function GOGraph. It contains
  #             all of the gene sets that will be in the graph.
  #   color.nodes: A character vector containing the gene sets that are of
  #                interest.
  #   interact: Indicates whether or not the graph should be
  #             interactive.
  #   legend.pos: Indicates what position the legend should be in.
  #   print.legend: Indicates whether or not a legend should be printed.
  #   use.col: Color to apply to nodes representing gene sets of interest
  #   bg.col: Color to use for border of all nodes, 
  #           names of all nodes NOT representing gene sets of interest, 
  #           and all edges
  nodes <- buildNodeList(GO.Graph)
  focusnode <- names(nodes) %in% color.nodes
  names(focusnode) <- names(nodes)
  nodefill <- ifelse(focusnode, use.col, "white")
  nAttrs <- list() 
  nAttrs$fillcolor <- nodefill
  nAttrs$label <- 1:length(names(nodes)) 
  names(nAttrs$label) <- names(nodes)
  #!# section added 07.11.14
  nAttrs$color <- rep(bg.col, length(nodes))
  names(nAttrs$color) <- names(nodes)
  fc <- ifelse(focusnode, 'black', bg.col)
  nAttrs$fontcolor <- fc
  names(nAttrs$fontcolor) <- names(nodes)
  edges <- buildEdgeList(GO.Graph)
  eNames <- names(edges)  
  eAttrs <- list()
  eAttrs$color <- rep(bg.col, length(eNames))
  names(eAttrs$color) <- eNames
  #!# section added 10.02.14
  eAttrs$arrowhead <- rep('none', length(eNames))
  names(eAttrs$arrowhead) <- eNames
  #!#
  pg <- plot(GO.Graph, nodeAttrs = nAttrs, edgeAttrs=eAttrs)
  x <- getNodeXY(pg)$x 
  y <- getNodeXY(pg)$y
  ordering <- sort.list(order(-y, x)) 
  nAttrs$label <- ordering
  names(nAttrs$label) <- names(nodes) 
  plot(GO.Graph, nodeAttrs = nAttrs, edgeAttrs=eAttrs)
  Terms <- sapply(lookUp(names(nodes)[sort.list(ordering)], "GO", "TERM"), Term) 
  names(Terms) <- NULL 
  legend <- data.frame(Terms)
  if(print.legend){
    print(legend)
  }
  if(interact){ 
    repeat {
      p <- locator(n = 1) 
      if (is.null(p)) break()
      pg <- plot(GO.Graph, nodeAttrs = nAttrs, edgeAttrs=eAttrs)
      x <- getNodeXY(pg)$x 
      y <- getNodeXY(pg)$y
      distance <- abs(p$x - x) + abs(p$y - y) 
      idx <- which.min(distance)
      legend(legend.pos, legend=c(nAttrs$label[idx], names(focusnode)[idx],
                                  Term(lookUp(names(focusnode)[idx], 
                                              "GO", "TERM")[[1]])), bg = "white")
    } 
  } 
}

########################################################

go2GeneSet <- function(name, object){
  # Creates a table showing the profile for a single gene set in each of the 
  # tested contrasts.
  #
  # Args:
  #   name: The name, or ID, of the desired gene set.
  #   object: A mvGST object containing a final results table.
  #
  # Returns:
  #   A matrix with possible profiles as the row names and contrasts as the
  #   column names. Ones in the appropriate cells showing which profile the
  #   gene set fit for each contrast and zeroes elsewhere.
  table <- object$results.table
  single.table <- matrix(0, nrow = nrow(table), ncol = ncol(table), dimnames = dimnames(table))
  observed <- profileCombine(object$ones.zeroes)
  profs <- lapply(observed, function(x) x[object$group.names == name, ])
  all.profs <- matrix(NA, nrow = nrow(table), ncol = ncol(observed[[1]]))
  nrt <- nrow(table)
  for (i in seq_len(nrt)){
    raw.profile <- dimnames(table)[[1]][i]
    temp <- substr(raw.profile, 3, nchar(raw.profile) - 1)
    temp2 <- as.integer(strsplit(temp, ",")[[1]])
    the.profile <- matrix(temp2, nrow = 1)
    all.profs[i, ] <- the.profile
  }
  ncs <- ncol(single.table)
  for (i in seq_len(ncs)){  
    if(is.null(dim(profs[[i]]))){   #!# Added this check in case a requested GO term wasn't included in gene sets tested
    row <- apply(t(profs[[i]] == t(all.profs)), MARGIN = 1, all)
    single.table[row, i] <- 1
	}
  }
  result <- list(results.table = single.table, ord.lev = object$ord.lev)
  class(result) <- "mvGST"
  return(single.table)
}

###########################################

fillInList <- function(group, term, offspring, list.groups){
  # Takes a gene sets in which genes are not associated
  # with offspring terms from the GO ontology and includes all genes
  # from offspring terms.
  #
  # Args:
  #   group: A character vector with the genes already associated
  #          with the set.
  #   term: The name of the gene set.
  #   offspring: A list showing the offspring sets of each gene set
  #   list.groups: A list showing all genes already associated with
  #                each gene set.
  if (is.na(offspring[[term]][1])){
    return(as.character(unlist(list.groups[term])))
  } else {
    new.group <- unique(unlist(c(group, list.groups[offspring[[term]]])))
    return(new.group)
  }
}

############################################

distributeWeight <- function(currentSets, children, weights, parentWeight = 1) 
{
  currentChildren <- children[currentSets]
  ready <- FALSE
  ii = 1
  while (!ready) {
    m <- length(currentSets)
    check <- weights
    weights[currentSets] <- weights[currentSets]+1/m*parentWeight
    if (length(weights)!=length(check) || class(weights)!='numeric') {
      warning("Weights changed length or class!!")
      flush.console()
      return(list(weights=weights, check=check, 
                  currentSets=currentSets, parentWeight=parentWeight, ii=ii))
    }
    if (ii > length(currentChildren)) {
      ii = 1
      currentSets <- unique(unlist(currentChildren))
      index <- currentSets %in% names(weights)
      if (sum(index) > 0) {
        currentChildren <- children[currentSets[index]]
        currentSets <- intersect(currentChildren[[ii]], names(weights))
        parentWeight <- weights[names(currentChildren[ii])]
        names(parentWeight) <- NULL
        ii = 2
      }
      else {
        ready <- TRUE
      }
    }
    else {
      currentSets <- intersect(currentChildren[[ii]], names(weights))
      parentWeight <- weights[names(currentChildren[ii])]  
      names(parentWeight) <- NULL
      ii = ii+1
    }
    
  }
  weights <- weights[!is.na(weights)]
  return(weights)
}

getCurrentChildren <- function(parents,sets) 
{
  tmp <- lapply(parents, function(pars) {
    if (sum(pars %in% sets)>0) TRUE
    else FALSE
  })
  tmp <- unlist(tmp)
  names(tmp[tmp])
}

makeCoherent <- function(currentSets,parents,weights,adjustedP)
{
  currentChildren <- getCurrentChildren(parents,currentSets)
  if (length(currentChildren)>0) {
    ii <- 1
    ready <- FALSE
    while (!ready) {
      currentChild <- currentChildren[ii]
      parentP <- min(adjustedP[parents[[currentChild]]])
      adjustedP[currentChild] <- max(parentP,adjustedP[currentChild])
      if (ii < length(currentChildren)) {
        ii = ii+1
      } 
      else {
        currentChildren <- getCurrentChildren(parents,currentChildren)
        if (length(currentChildren) < 1) 
          ready <- TRUE
        else
          ii <- 1
      }
    }
  }
  return(adjustedP)
}


getAncestorsAndOffspring <- function(GOid,ontology = c('MF','CC','BP')) {
  ancestors <- lapply(as.list(ontology), function(ont) {
    ext <- paste(ont, "ANCESTOR", sep = "")
    GOOBJ <- eval(as.name(paste("GO", ont, "ANCESTOR", 
                                sep = "")))
    ontid <- intersect(keys(GOOBJ), GOid)
    if (length(ontid) > 0) {
      tmp <- lookUp(ontid, "GO", ext)
      choose <- which(names(tmp)%in%GOid)
      tmp[choose]
      tmp <- lapply(tmp,function(t) {
        intersect(setdiff(t,"all"),
                  GOid)})
    } else list()
  })
  ancestors <- do.call(c, ancestors)
  offspring <- lapply(as.list(ontology), function(ont) {
    ext <- paste(ont, "OFFSPRING", sep = "")
    GOOBJ <- eval(as.name(paste("GO", ont, "OFFSPRING", 
                                sep = "")))
    ontid <- intersect(keys(GOOBJ), GOid)
    if (length(ontid) > 0) {
      tmp <- lookUp(ontid, "GO", ext) 
      choose <- which(names(tmp)%in%GOid)
      tmp[choose]
      tmp <- lapply(tmp,function(t) intersect(t,GOid))
    } else list()
  })
  offspring <- do.call(c, offspring)
  offspring <- sapply(offspring, function(os) 
    if (all(is.na(os))) character(0)  else os)
  
  return(list( GOid = GOid, ancestors = ancestors, 
               offspring = offspring ))
}


turnListAround <- function (aList) #from globaltest package
{
  newlist <- new.env(hash = TRUE)
  objs <- names(aList)
  if (is.null(objs)) 
    objs <- seq_along(alist)
  for (i in objs) {
    for (j in aList[[i]]) {
      newlist[[j]] <- c(newlist[[j]], i)
    }
  }
  as.list(newlist)
}