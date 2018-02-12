profileTable <-
function(pvals, gene.names = NULL, contrasts = NULL, list.groups = NULL, sig.level = .05, 
                         gene.ID, organism, affy.chip, ontology = "BP", method = 2,
         minsize = 1, maxsize = Inf, mult.adj = "BY"){
  # Takes a matrix of p-values and returns a table with the number of gene sets 
  # that fall into each possible profile
  #
  # Args:
  #   gene.names: A character vector containing the gene names that correspond
  #               to the rows of the matrix of p-values
  #   contrasts: A character vector containing the contrasts that correspond to
  #              each column in the matirx of p-values. Must be in 1 of 2 formats:
  #              Var1.Var2 or Var1
  #   pvals: A matrix containing the p-values corresponding to the various genes
  #          and contrasts
  #   list.groups: An optional list containing user-defined gene sets
  #   sig.level: The alpha level that should be used. Default is .05.
  #   gene.ID: Gene naming system used for the gene names. Used to generate list
  #            of gene sets mapping genes to Gene Ontology sets
  #   organism: The organism that the genes come from. Used to generate list
  #            of gene sets mapping genes to Gene Ontology sets
  #   affy.chip: The type of affy.chip used, if gene.ID == "affy".
  #   ontology: The ontolgy that should be used for gene sets: BP, MF, or CC 
  #   method: The method for handling gene name tranlation issues
  #   minsize: The minimum size gene set that will be included in the list
  #   maxsize: The maximum size gene set that will be included in the list
  #   mult.adj: "BY" for Benjamini-Yekutieli adjustment. "SFL" for short focus level
  #             adjustment. "none" for no adjustment.  Default is "BY".
  #
  # Returns:
  #   The primary return is a table containing the number of gene sets that fall
  #   into each of the possible significance profiles
  if (is.null(contrasts)){
    contrasts <- colnames(pvals)
  }
  if (is.null(gene.names)){
    gene.names <- rownames(pvals)
  }
  if (!is.matrix(pvals)){
    stop("pvals must be a matrix")
  }
  if(!is.element(mult.adj,c("BY","SFL"))){
    stop("mult.adj must be one of 'BY' or 'SFL'")
  }
  # adding ".BP" to contrasts if in the form Var1
  vars <- length(unlist(strsplit(contrasts[1], "[.]")))
  if (vars == 1){
    contrasts <- paste(contrasts, ontology, sep = ".")
  }
  # accounting for the possibilities of 1's and 0's in the p-values
  f.one <- 1-.Machine$double.eps
  f.zero <- .Machine$double.eps
  pvals <- ifelse(pvals == 1, f.one, pvals)
  pvals <- ifelse(pvals == 0, f.zero, pvals)
  
  if (is.null(list.groups)){
    if (any(gene.ID == c("affy", "genbank", "alias", "ensembl", "entrez", 
                         "symbol", "genename", "unigene")) != TRUE){  
	   # Converts gene names to entrez database that can be mapped directly to GO.
        old.names <- character()
        new.names <- character()
        # gene names are converted 1000 at a time because doing more than that
        # may cause an error (vector too large)
		for(i in seq.int(from=0, to=floor(length(gene.names)/1000))){
        #for (i in 0:floor(length(gene.names) / 1000)){
          low <- i * 1000 + 1
          high <- (i + 1) * 1000
          if (high > length(gene.names)){ 
            high <- length(gene.names)
          }
          gene.names <- toupper(gene.names)
          converted1 <- gconvert(gene.names[low:high], target = "ENTREZGENE_ACC", organism=organism)
          old.names <- c(old.names, as.character(converted1$alias))
          new.names <- c(new.names, as.character(converted1$target))
        }

      # trimming all-numeric ID's down to just the numbers
      #!#if(length(grep("^[[:digit:]]+$", str_trim(gene.names[1]))) == 1){ 
      #!#  all.numeric <- TRUE
      #!#}else{all.numeric <- FALSE} #!# added length() and else{} to this
      if(grep("^[[:digit:]]+$", str_trim(gene.names[1])) == 1){ 
          all.numeric <- TRUE
        }


      new.genes <- cbind(old.names, new.names)
      new.genes[, 2] <- gsub("ENTREZGENE_ACC:", "", new.genes[, 2])
      if (all.numeric){
        new.genes[, 1] <- substr(new.genes[, 1], 
                                 regexpr("[[:digit:]]+$", new.genes[, 1]),
                                 nchar(new.genes[, 1]))
      }
      # eliminating duplicate rows (i.e., both the new and old gene names
      #                             are the same)
      duplicate <- rep(FALSE, length(new.genes[, 1]))
      for (i in seq_along(new.genes[, 1])[-1]){
        duplicate[i] <- ifelse(all(new.genes[i, ] == new.genes[i-1, ]), TRUE, FALSE)
      }
      new.genes <- new.genes[!duplicate, ]
      
      converted <- geneNameConvertRows(pvals, gene.names, new.genes, method)
      pvals <- converted$p.mat
      gene.names <- converted$genes
      gene.ID <- "entrez"
    } 
    # generates list of gene sets
    list.groups1 <- generateGeneSets(ontology=ontology, species = organism, 
                                    ID = gene.ID, affy.chip)
    offspring <- get("as.list", pos = "package:AnnotationDbi")(get(paste("GO", ontology, "OFFSPRING", sep = "")))
    list.groups <- sapply(1:length(offspring), 
                          function(x) fillInList(list.groups1[[names(offspring[x])]], 
                                                 names(offspring)[x],
                                                 offspring, list.groups1))
    names(list.groups) <- names(offspring)
    size <- sapply(list.groups, length)
    #list.groups <- list.groups[keepers] #!# commented out this line; keepers apparently an artifact from method4
    list.groups <- list.groups[(size>=minsize)&(size<=maxsize)]  #!# added in version 0.99.0 to allow gene set size restrictions
  }
  # Creates MVGST object
  pmat <- mvGSTObject(gene.names, contrasts, pvals, list.groups)
  # Converts p-values matrix for genes to a p-value matrix for gene sets
  grouped.pmat <- separate(pmat, list.groups)
  
  # Perform multiple hypothesis test adjustment for multiple hypothesis testing
  if (mult.adj == "SFL"){
    adjusted <- grouped.pmat
    adjusted$adjusted.group.pvals <- grouped.pmat$grouped.raw
    adjusted$adjusted.group.pvals[] <- NA
	ncgp <- ncol(grouped.pmat$grouped.raw)
    for (i in seq_len(ncgp)){
      pvalues <- grouped.pmat$grouped.raw[, i]
      two.sided <- convertPvalues(pvalues, two.sided = FALSE)
      names(two.sided) <- grouped.pmat$group.names
      new.pvalues <- p.adjust.SFL(two.sided, sig.level=sig.level)
      relative <- ifelse(pvalues < .5, 1, -1)
      one.combined <- convertPvalues(new.pvalues, relative)
      adjusted$adjusted.group.pvals[, i] <- one.combined
    }
  } else if(mult.adj =="BY"){
    adjusted <- oneSideBYAdjust(grouped.pmat)
      } 
	  
  # Creates the final output
  one.zero <- changeTO10(adjusted, sig.level=sig.level)
  almost.final <- finalResults(one.zero)
  near.final <- cut(almost.final)
  final <- mvSort(near.final)
  return(final)
}
