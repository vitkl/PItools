##' Filter molecular interactions by Pubmed ID
##' @name subsetMITABbyPMIDs
##' @author Vitalii Kleshchevnikov
##' @description filter molecular interaction data by Pubmed ID of publications
##' @param MITABdata object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by  \code{\link[PItools]{cleanMITAB}})
##' @param PMIDs character or character vector, return only interactions reported by these Pubmed IDs
##' @param inverse_filter logical, inverse filtering criteria
##' @return return only interactions reported by \code{PMIDs}
##' @import data.table
##' @seealso \code{\link[PItools]{subsetMITABbyMethod}}), \code{\link[PItools]{subsetMITABbyID}})
##' @export subsetMITABbyPMIDs
##' @examples
##' {
##' # Download all human interaction data
##' full = fullInteractome(taxid = "9606", database = "IntActFTP",
##'          clean = TRUE, protein_only = TRUE)
##'
##' # subset unpublished Vidal group data
##' Vidal = subsetMITABbyPMIDs(MITABdata = full,
##'                PMIDs = "unassigned1304")
##'
##' # subset both published and unpublished Vidal group data
##' Vidal_all = subsetMITABbyPMIDs(MITABdata = full,
##'                PMIDs = c("25416956", "unassigned1304"))
##'
##' # subset Mattias Mann 2015 paper data
##' Mann = subsetMITABbyPMIDs(MITABdata = full,
##'                PMIDs = "26496610")
##' }
subsetMITABbyPMIDs = function(MITABdata, PMIDs = NULL, inverse_filter = F) {
  if(!grepl("clean_MItab",class(MITABdata))) stop("MITABdata is not of class clean_MItab27 or related clean_MItab class")
  if(!is.null(PMIDs) | is.character(PMIDs) | is.integer(PMIDs) | is.numeric(PMIDs)) {
    NULL
  } else {
    stop("PMIDs not supplied or not a vector")
  }

  if(inverse_filter){
    MITABdata$data = MITABdata$data[!Publication_Identifiers %in% PMIDs,]
  } else MITABdata$data = MITABdata$data[Publication_Identifiers %in% PMIDs,]

  MITABdata$subsetByPMIDsDetails = MITABdata$data[, table(Publication_Identifiers)]
  warn = MITABdata$subsetByPMIDsDetails == 0
  if(sum(warn) > 0) warning(paste0("no data for the following PMIDs: ", paste0(names(MITABdata$subsetByPMIDsDetails[warn]),collapse = ", ")))

  return(MITABdata)
}

printSubsetByPMIDsDetails = function(data){
  cat(paste0("\n` The data comes from the following publications: `\n"))
  print(data$subsetByPMIDsDetails)
}
