pickOut <-
function(mvgst, row, col = 1){
  # Return the gene sets that are included in a specific cell of the final
  # results table of a mvGST object.
  #
  # Args:
  #   mvgst: A mvGST object with a final results table.
  #   row: The row of the desired cell.
  #   column: The column of the desired cell.
  #
  # Returns:
  #   A character vector containing the gene sets in the desired cell.
  table <- mvgst$results.table
  raw.profile <- dimnames(table)[[1]][row]
  #!# added following section 07.10.14
  splt.raw.profile <- unlist(strsplit(raw.profile, ''))
  if(splt.raw.profile[1] == "c")  # i.e., if(profile is multivariate)
    {
      temp <- substr(raw.profile, 3, nchar(raw.profile) - 1)
      temp2 <- as.integer(strsplit(temp, ",")[[1]])
      the.profile <- matrix(temp2, nrow = 1)
     }else{
           the.profile <- matrix(raw.profile, nrow=1)  
           }
  #!# end of 07.10.14 edit
  observed <- profileCombine(mvgst$ones.zeroes)[[col]] 
  t <- apply(observed, MARGIN = 1, FUN = function(x) all(x == the.profile))
  GO.ID <- mvgst$group.names[!is.na(t) & t]
  GO.Description <- unname(unlist(getGOTerm(GO.ID)))
  index <- mvgst$group.names %in% GO.ID
  adj.pvals <- mvgst$adjusted.group.pvals[index,]
  rownames(adj.pvals) <- NULL
  output <- as.data.frame(cbind(GO.ID, GO.Description, adj.pvals))
  return(output)
}
