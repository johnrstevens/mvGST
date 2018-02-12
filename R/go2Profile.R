go2Profile <-
function(names, object){
  # Performs the same operation as go2GeneSet, but for multiple gene sets
  #
  # Args:
  #   names: A character vector with the names, ID's, of the gene sets
  #          of interest.
  #   object: A mvGST object with a final results table.
  #
  # Returns:
  #   A list of matrices. Each matrix has possible profiles as the 
  #   row names and contrasts as the column names. Ones in the appropriate 
  #   cells showing which profile the gene set fit for each contrast and 
  #   zeroes elsewhere.
  result <- lapply(names, go2GeneSet, object)
  list.results <- list()
  temp <- list(results.table = NULL, ord.lev = NULL)
  temp$ord.lev <- object$ord.lev
  for (i in seq_along(result)){
    temp$results.table <- result[[i]] 
    class(temp) <- "mvGST"
    temp1 <- cut(temp)
    list.results[[i]] <- temp1
  }
  names(list.results) <- names
  return(list.results)
}
