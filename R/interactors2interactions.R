##' Retrieve interactions for a list of interactors (given data table, not using webservice)
##' @name interactors2interactions
##' @author Vitalii Kleshchevnikov
##' @description subset molecular interaction data (cleaned MITAB format in a data.table object) with a list of interactors
##' @param MITABdata data.table with the pair_id column - unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs
##' @param interactorsList filter \code{MITABdata} using this list of interactors
##' @return data.table, a subset of MITABdata that contains interactions for interactorsList molecules
##' @import data.table
##' @export interactors2interactions
interactors2interactions = function(MITABdata, interactorsList){
  interactionsubset = copy(MITABdata)
  interactionsubset[, c("ida", "idb") := tstrsplit(pair_id, "\\|")]
  interactionsubset = interactionsubset[ida %in% interactorsList | idb %in% interactorsList, ]
  interactionsubset[, c("ida", "idb") := NULL]
  return(interactionsubset)
}
