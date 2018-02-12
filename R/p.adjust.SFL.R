p.adjust.SFL <- function (rawp, ontology = c("BP", "CC", "MF"), focus = "rn", 
    ancestors, offspring, trace = FALSE, recycle = TRUE, sig.level = 0.05) 
{
    if (missing(rawp)) 
        stop("argument \"rawp\" is missing with no default")
    if (is.null(names(rawp))) 
        stop("argument \"rawp\" must be named")
    else GOid <- names(rawp)
    myGOTERM <- lookUp(GOid, "GO", "TERM")
    if (length(ontology) > 1) 
        ontology <- "BP"
    choose <- sapply(myGOTERM, Ontology) %in% ontology
    if (sum(choose) < length(rawp)) {
        warning("Some values of \"rawp\" were not in the correct ontology and were removed.")
        if (trace) 
            which(!GOid %in% GOid[choose])
        GOid <- GOid[choose]
        rawp <- rawp[GOid]
    }
    if (missing(ancestors) && missing(offspring)) {
        if (trace) {
            message("Getting ancestors and offspring: ", date())
            flush.console()
        }
        ancestors <- getAncestorsAndOffspring(GOid, ontology)
        offspring <- ancestors$offspring
        ancestors <- ancestors$ancestors
    } 
    if (missing(ancestors)) 
        ancestors <- turnListAround(offspring)
    if (missing(offspring)) 
        offspring <- turnListAround(ancestors)
	## begin section added 8/26/14 to have better error catching and restrict GO terms to names in rawp
    ancestors <- ancestors[GOid]
    offspring <- offspring[GOid]
    ancestors <- lapply(ancestors, function(anc){
        intersect(anc,GOid)
    })
    offspring <- lapply(offspring, function(off){
        intersect(off,GOid)
    })
    if (length(focus) > 1)
        focus <- intersect(focus,GOid)
	## end section added 8/26/14
    else if (focus == "rn") 
        focus <- names(which.max(lapply(offspring, length)))
    else
        stop("focuslevel must be a character vector of GO IDs or \"rn\"") 
    message("focus (length ", length(focus),") = ", focus)
    if (trace) {
        message("Adjusting raw p-values: ", date())
        flush.console()
    }
    rawp[is.na(rawp)] <- 1
    if (trace) {
        message("Identifying child relations: ", date())
        flush.console()
    }
    children <- lapply(offspring, function(offs) {
        setdiff(offs, unique(unlist(offspring[offs])))
    })
    currentSets <- focus
    weights <- rep(0, length(GOid))
    names(weights) <- GOid
    if (trace) {
        message("Distributing weights: ", date())
        flush.console()
    }
    weights <- distributeWeight(currentSets, children, weights)
    adjustedP <- rawp/weights
    index <- adjustedP > 1
    adjustedP[index] <- 1
    if (trace) {
        message("Identifying parent relations: ", date())
        flush.console()
    }
    parents <- lapply(ancestors, function(anc) {
        setdiff(anc, unique(unlist(ancestors[anc])))
    })
    if (trace) {
        message("Ensuring coherence of p-values: ", date())
        flush.console()
    }
    adjustedP <- makeCoherent(focus, parents, weights, adjustedP)
    if (recycle) {
        if (trace) {
            message("Recycling thresholds: ", date())
            flush.console()
        }
        ready <- FALSE
        setsFocusDown <- unique(unlist(offspring[focus]))
        parentsFocusDown <- parents[setsFocusDown]
        if (!length(parentsFocusDown) > 0) 
            parentsFocusDown <- NULL
        setsFocusDown <- union(focus, setsFocusDown)
        barren <- sapply(children, length) == 0
        EndNodes <- names(barren[barren])
        sigEndNodesUsed <- logical(length(EndNodes))
        sigEndNodes <- adjustedP[EndNodes] <= sig.level
        while (!ready) {
            if (all(!sigEndNodes[!sigEndNodesUsed])) {
                ready <- TRUE
            }
            else {
                freeWeight <- sum(weights[EndNodes][sigEndNodes])
                nonRejected <- setsFocusDown[adjustedP[setsFocusDown] > 
                  sig.level]
                tmp <- sapply(nonRejected, function(set) {
                  index <- parents[[set]] %in% nonRejected
                  if (!all(index)) 
                    names(parents[set])
                })
                topNodes <- unique(unlist(tmp))
                weights[setdiff(nonRejected, topNodes)] <- 0
                tmp <- distributeWeight(topNodes, children, weights[nonRejected], 
                  freeWeight)
                weights[nonRejected] <- tmp
                adjustedP[nonRejected] <- rawp[nonRejected]/weights[nonRejected]
                index <- adjustedP > 1
                adjustedP[index] <- 1
                adjustedP <- makeCoherent(focus, parents, weights, 
                  adjustedP)
                sigEndNodesUsed <- sigEndNodes
                sigEndNodes <- adjustedP[EndNodes] <= sig.level
            }
        }
    }
    rootfocus <- names(which.max(lapply(offspring, length)))
    if (length(setdiff(rootfocus, focus)) > 0) {
        focusAncestors <- unique(unlist(ancestors[focus]))
        adjustedP[focusAncestors] <- sapply(focusAncestors, function(anc) min(adjustedP[offspring[[anc]] %in% 
            focus]))
    }
    adjustedP
}
