graphCell <-
function(object, row, col = 1, set = NULL, ontology = "BP", interact = TRUE, legend.pos = "bottomleft",
                      print.legend = TRUE, use.col="red", bg.col = "grey80"){
  # Graphs the gene sets from a given cell of the final results table
  # of a mvGST object.
  #
  # Args:
  #   object: A mvGST object with a final results table.
  #   row: The row of the desired cell.
  #   col: The column of the desired cell.
  #   set: optional, data frame with the first column containing the gene sets
  #         that should be in the GO graph.
  #   ontology: The ontology, within Gene Ontology, that should be used ("BP", "MF", "CC").
  #   interact: Indicates whether or not the graph should be interactive.
  #   legend.pos: If interactive, indicates the desired position of the legend.
  #   print.legend: Indicates if the legend should also be printed separately
  #   use.col: Color to apply to nodes representing gene sets of interest
  #   bg.col: Color to use for border of all nodes, 
  #           names of all nodes NOT representing gene sets of interest, 
  #           and all edges
  if (is.null(set)){
    set <- pickOut(object, row, col)
  } 
  set <- as.character(set[, 1])
  g.sub <- switch(ontology, 
                   BP = GOGraph(set,GOBPPARENTS),
                   MF = GOGraph(set,GOMFPARENTS),
                   CC = GOGraph(set,GOCCPARENTS),
				   stop("ontology must be one of 'BP','MF', or 'CC'")
				  )
  g.sub <- removeNode("all",g.sub)
  interactiveGraph(g.sub, set, interact, legend.pos, print.legend, use.col, bg.col)
  invisible(NULL)
}
