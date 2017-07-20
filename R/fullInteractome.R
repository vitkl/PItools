##' Retrieve interactome (proteins only or not) of a given taxonomic species from a particular database
##' @name fullInteractome
##' @author Vitalii Kleshchevnikov
##' @description Retrieve interactome (proteins only or not) of a given taxonomic species from a particular database. Interactome can be additionally cleaned to include only specific information: \code{\link{cleanMITAB}}
##' @details \code{taxid} is used to query specified database using PSICQUIC client, only interactions in which both participants belong the \code{taxid} are retured (\code{"taxidA:9606 AND taxidB:9606"}, not \code{"species:9606"}). Details: \code{\link{queryPSICQUIC}}
##' @details \code{fullInteractome} can be used to retrive interactome data using PSICQUIC service using \code{\link{queryPSICQUIC}}, clean and select specific columns using \code{\link{cleanMITAB}} and filter resulting dataset for protein-protein interaction only. This is the default option.
##' @details Alternatively, \code{fullInteractome} can only retrive interactome data using PSICQUIC service without cleaning of filtering.
##' @details Another option is to supply \code{MITABdata} to be cleaned and filtered
##' @param MITABdata data.table containing molecular interaction data as returned by \code{\link{queryPSICQUICrlib}}, default in NULL
##' @param taxid character (1L), taxonomy id of the species which interaction participants should belong to, default is "9606" (which is human)
##' @param database character (1L), argument for \code{\link{queryPSICQUIC}}, default is "imex"
##' @param format character (1L), argument for \code{\link{queryPSICQUIC}}, default is "tab25"
##' @param clean logical (1L), if TRUE extract specific information using \code{\link{cleanMITAB}}, default is TRUE
##' @param protein_only logical (1L), if TRUE the interaction participants are restricted to proteins (exclude other types of molecules such as RNA or small molecules), default is TRUE
##' @param directory directory where to store the data, if NULL the data is stored in <R-package-library>/MItools/data
##' @return data.table containing molecular interaction data in either of these two formats:
##' @return if \code{clean} is TRUE: contains columns as described in \code{\link{cleanMITAB}};
##' @return if \code{clean} is FALSE: contains a standard set of columns for MITAB2.5 or MITAB2.7 depending on \code{format};
##' @import data.table
##' @export fullInteractome
##' @examples
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.5 format, cleans and select specific columns
##' full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE)
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.5 format, cleans and select specific columns; save it to the specific directory inside working directory
##' full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE, directory = "./data/")
fullInteractome = function(MITABdata = NULL, taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE, directory = NULL){
  # if the interaction data for species taxid and from database is not saved in the library - queryPSICQUIC for interaction data for taxid interactions in the database and in MITAB2.5 format, save results to the library
  if(is.null(MITABdata)){
    full_interactome = queryPSICQUICrlib(query = paste0("taxidA:",taxid," AND ", "taxidB:",taxid),
                                         format = format,
                                         database = database,
                                         directory = directory)
  }
  if(!is.null(MITABdata)) full_interactome = copy(MITABdata)

  # clean this data to make it more useble if clean is TRUE
  if(clean){
    full_interactome_clean = cleanMITAB(full_interactome)
    # filter out non-proteins
    if(protein_only){
      # removing interactions if at least one interactor has non-uniprot id
      full_interactome_clean = full_interactome_clean[interactor_IDs_databases_A == "uniprotkb" & interactor_IDs_databases_B == "uniprotkb",]
    }
    return(full_interactome_clean)
  }
  # if clean is FALSE return the interactome date in the raw MITAB format
  if(!clean){
    return(full_interactome)
  }

}
