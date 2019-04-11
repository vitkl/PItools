##' Remove self loops (self-interactions)
##' @name removeSelfLoops
##' @author Vitalii Kleshchevnikov
##' @description subset molecular interaction data with a list of interactors
##' @param MITABdata molecular interaction data, object of class "clean_MItab" (or any other object which names starts from "clean_MItab")
##' @return molecular interaction data with the self interactions removed
##' @import data.table
##' @seealso \code{\link[PItools]{subsetMITABbyMethod}}), \code{\link[PItools]{subsetMITABbyPMIDs}})
##' @export removeSelfLoops
##'
removeSelfLoops = function(MITABdata){
  if(!grepl("clean_MItab",class(MITABdata))) stop("MITABdata is not of class clean_MItab27 or related clean_MItab class")
  MITABdata$data = MITABdata$data[!IDs_interactor_A == IDs_interactor_B,]
  MITABdata
}
