##' Retrieve interactions for a list of interactors (given data table, not using webservice)
##' @name subsetMITABbyID
##' @author Vitalii Kleshchevnikov
##' @description subset molecular interaction data with a list of interactors
##' @param MITABdata object of class "clean_MItab25", "clean_MItab27", "clean_MItab25_fullInteractome", "clean_MItab27_fullInteractome"
##' @param ID_seed filter \code{MITABdata} using this list of interactors
##' @param within_seed logical, should \code{subsetMITABbyID} return only interactions between molecules in \code{ID_seed} (TRUE) or any interactions molecules in \code{ID_seed} have in \code{MITABdata}
##' @return object of class `class of MITABdata`_subset (for example, clean_MItab27_subset, basically list), a subset of \code{MITABdata} that contains interactions for \code{ID_seed} molecules, and additionally includes the ID_seed and within_seed
##' @import data.table
##' @seealso \code{\link[MItools]{subsetMITABbyMethod}}), \code{\link[MItools]{subsetMITABbyPMIDs}})
##' @export subsetMITABbyID
##' @export print.clean_MItab25_subset
##' @export print.clean_MItab27_subset
##' @export print.clean_MItab25_fullInteractome_subset
##' @export print.clean_MItab27_fullInteractome_subset
##' @export print.clean_MItab25_interSpeciesInteractome_subset
##' @export print.clean_MItab27_interSpeciesInteractome_subset
subsetMITABbyID = function(MITABdata, ID_seed, within_seed = F, only_seed2nonseed = F){
  if(!grepl("clean_MItab",class(MITABdata))) stop("MITABdata is not of class clean_MItab27 or related clean_MItab class")

  if(within_seed) {
    if(only_seed2nonseed) stop("you can select interactions both only within seed proteins (within_seed = T) AND only between seed and non-seed proteins (only_seed2nonseed = T)")
    MITABdata$data = MITABdata$data[IDs_interactor_A %in% ID_seed & IDs_interactor_B %in% ID_seed, ]
    }
  if(!within_seed) MITABdata$data = MITABdata$data[IDs_interactor_A %in% ID_seed | IDs_interactor_B %in% ID_seed, ]

  if(only_seed2nonseed){
    MITABdata$data = MITABdata$data[(IDs_interactor_A %in% ID_seed & !(IDs_interactor_B %in% ID_seed)) |
                                      (IDs_interactor_B %in% ID_seed & !(IDs_interactor_A %in% ID_seed)), ]
    MITABdata$data[IDs_interactor_A %in% ID_seed, IDs_A_order := IDs_interactor_A]
    MITABdata$data[IDs_interactor_A %in% ID_seed, IDs_B_order := IDs_interactor_B]
    MITABdata$data[!(IDs_interactor_A %in% ID_seed), IDs_A_order := IDs_interactor_B]
    MITABdata$data[!(IDs_interactor_A %in% ID_seed), IDs_B_order := IDs_interactor_A]
    if(grepl("clean_MItab25",class(MITABdata))) MITABdata$data = reorderMITAB25(MITABdata$data)
    if(grepl("clean_MItab27",class(MITABdata))) MITABdata$data = reorderMITAB27(MITABdata$data)
  }

  MITABdata$ID_seed = ID_seed
  MITABdata$within_seed = within_seed
  MITABdata$only_seed2nonseed = only_seed2nonseed

  if(!grepl("_subset",class(MITABdata))) class(MITABdata) = paste0(class(MITABdata), "_subset")
  return(MITABdata)
}

#print methods
print.clean_MItab25_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_subset, which is a subset of the data for query, file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab27_subset, which is a subset of the data for query, file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}

print.clean_MItab25_fullInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_fullInteractome_subset, which is a subset of the full interactome for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_fullInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab27_fullInteractome_subset, which is a subset of the full interactome for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}

print.clean_MItab25_interSpeciesInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_interSpeciesInteractome_subset, which is a subset of the clean_MItab25_interSpeciesInteractome that contains interactions between molecules (proteins, RNA) of taxid1: ", data$taxid1, " and taxid2: ", data$taxid2," , proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions, seed proteins are in the IDs_interactor_A)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_interSpeciesInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_interSpeciesInteractome_subset, which is a subset of the clean_MItab25_interSpeciesInteractome that contains interactions between molecules (proteins, RNA) of taxid1: ", data$taxid1, " and taxid2: ", data$taxid2," , proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, \ndoes it include ONLY interactions between these molecules (within_seed)?: ", data$within_seed), " \n does it includes only interactions between seed and non-seed proteins (excluding within seed interactions, seed proteins are in the IDs_interactor_A)?: ",data$only_seed2nonseed,"`\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}
