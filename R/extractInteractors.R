##' Subset MITAB to retrieve interactors, the number of interactions ot interactors
##' @rdname extractInteractors
##' @name extractInteractors
##' @author Vitalii Kleshchevnikov
##' @description subset molecular interaction data (cleaned MITAB format in a data.table object) with a list of interactors
##' @param cleanMITAB object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by  \code{\link[MItools]{cleanMITAB}})
##' @param taxid filter \code{cleanMITAB} using this list of taxonomy IDs
##' @param inverse_filter logical, inverse filtering criteria
##' @return character vector of interactor IDs or an integer count of interactions or interactors
##' @import data.table
##' @export extractInteractors
##' @export NuniqueInteractions
##' @export NuniqueInteractors
extractInteractors = function(cleanMITAB, taxid = NULL, inverse_filter = F) {
if(!grepl("clean_MItab",class(cleanMITAB))) stop("cleanMITAB is not of class clean_MItab27 or related clean_MItab class")
if(is.null(taxid)){
  res = cleanMITAB$data[, unique(c(IDs_interactor_A, IDs_interactor_B))]
} else {
  if(inverse_filter){
    res = cleanMITAB$data[!Taxid_interactor_A %in% taxid, unique(IDs_interactor_A)]
    res = c(res, cleanMITAB$data[!Taxid_interactor_B %in% taxid, unique(IDs_interactor_B)])
    res = unique(res)
  } else {
    res = cleanMITAB$data[Taxid_interactor_A %in% taxid, unique(IDs_interactor_A)]
    res = c(res, cleanMITAB$data[Taxid_interactor_B %in% taxid, unique(IDs_interactor_B)])
    res = unique(res)
  }
}
}

##' @rdname extractInteractors
##' @name NuniqueInteractions
##' @import data.table
##' @export NuniqueInteractions
NuniqueInteractions = function(cleanMITAB){
  if(!grepl("clean_MItab",class(cleanMITAB))) stop("cleanMITAB is not of class clean_MItab27 or related clean_MItab class")
  length(unique(cleanMITAB$data$pair_id))
}
##' @rdname extractInteractors
##' @name NuniqueInteractors
##' @import data.table
##' @export NuniqueInteractors
NuniqueInteractors = function(cleanMITAB, taxid = NULL, inverse_filter = NULL){
  if(!grepl("clean_MItab",class(cleanMITAB))) stop("cleanMITAB is not of class clean_MItab27 or related clean_MItab class")

  if(is.null(taxid)){
    return(length(cleanMITAB$data[,unique(c(IDs_interactor_A, IDs_interactor_B))]))
  } else {
    proteins = extractInteractors(cleanMITAB, taxid = taxid, inverse_filter = inverse_filter)
    return(length(proteins))
  }
}
