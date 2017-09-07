##' Retrieve interactions for a list of interactors (given data table, not using webservice)
##' @name subsetMITABbyID
##' @author Vitalii Kleshchevnikov
##' @description subset molecular interaction data (cleaned MITAB format in a data.table object) with a list of interactors
##' @param MITABdata object of class "clean_MItab25", "clean_MItab27", "clean_MItab25_fullInteractome", "clean_MItab27_fullInteractome"
##' @param ID_seed filter \code{MITABdata} using this list of interactors
##' @param within_seed logical, should \code{subsetMITABbyID} return only interactions between molecules in \code{ID_seed} (TRUE) or any interactions molecules in \code{ID_seed} have in \code{MITABdata}
##' @return object of class `class of MITABdata`_subset (for example, clean_MItab27_subset, basically list), a subset of \code{MITABdata} that contains interactions for \code{ID_seed} molecules, and additionally includes the ID_seed and within_seed
##' @import data.table
##' @export subsetMITABbyID
##' @export print.clean_MItab25_subset
##' @export print.clean_MItab27_subset
##' @export print.clean_MItab25_fullInteractome_subset
##' @export print.clean_MItab27_fullInteractome_subset
subsetMITABbyID = function(MITABdata, ID_seed, within_seed = F){
  valid_class = c("clean_MItab25", "clean_MItab27", "clean_MItab25_fullInteractome", "clean_MItab27_fullInteractome")
  if(!(class(MITABdata) %in% valid_class)) stop("subsetMITABbyID works only on objects of class: ", paste(valid_class, collapse = ", "))

  if(within_seed) MITABdata$data = MITABdata$data[IDs_interactor_A %in% ID_seed & IDs_interactor_B %in% ID_seed, ]
  if(!within_seed) MITABdata$data = MITABdata$data[IDs_interactor_A %in% ID_seed | IDs_interactor_B %in% ID_seed, ]

  MITABdata$ID_seed = ID_seed
  MITABdata$within_seed = within_seed

  if(!grepl("_subset",class(MITABdata))) class(MITABdata) = paste0(class(MITABdata), "_subset")
  return(MITABdata)
}

#print methods
print.clean_MItab25_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_subset, which is a subset of the data for query, file, format, databases, date: `\n"))
  print(data$metadata)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, does it include ONLY interactions between these molecules? ", data$within_seed), " `\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab27_subset, which is a subset of the data for query, file, format, databases, date: `\n"))
  print(data$metadata)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, does it include ONLY interactions between these molecules? ", data$within_seed, " `\n"))
  cat("\n` view of the $data: `\n")
  print(data$data)
}

print.clean_MItab25_fullInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab25_fullInteractome_subset, which is a subset of the full interactome for taxid: ", data$taxid, ", proteins only: ", protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, does it include ONLY interactions between these molecules? ", data$within_seed, " `\n"))
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_fullInteractome_subset = function(data){
  cat(paste0("\n` Object of class clean_MItab27_fullInteractome_subset, which is a subset of the full interactome for taxid: ", data$taxid, ", proteins only: ", protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat(paste0("\n` this subset contains interactions of a set of molecules (seed): ", length(data$ID_seed), " total, does it include ONLY interactions between these molecules? ", data$within_seed)," `\n")
  cat("\n` view of the $data: `\n")
  print(data$data)
}