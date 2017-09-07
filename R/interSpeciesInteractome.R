##' Retrieve interactome (proteins only or not) of a given pair of taxonomic species from a particular database
##' @name interSpeciesInteractome
##' @author Vitalii Kleshchevnikov
##' @description Retrieve interactome (proteins only or not) of a given pairs of taxonomic species from a particular database. Interactome can be additionally cleaned and includes only specific information: \code{\link{cleanMITAB}}
##' @details \code{taxid1} and \code{taxid2} is used to query specified database using PSICQUIC client, only interacting pairs between \code{taxid1} and \code{taxid2} are retured (no interactions within the same species, \code{"(taxidA:taxid1 AND taxidB:taxid2) OR (taxidA:taxid2 AND taxidB:taxid1)"}).
##' @details \code{interSpeciesInteractome} can be used to retrive interactome data using PSICQUIC service using \code{\link{queryPSICQUIC}}, clean and select specific columns using \code{\link{cleanMITAB}} and filter resulting dataset for protein-protein interaction only. This is the default option.
##' @details Alternatively, \code{interSpeciesInteractome} can only retrive interactome data using PSICQUIC service without cleaning of filtering.
##' @details Another option is to supply \code{MITABdata} to be cleaned and filtered
##' @details Finally, you can avoid using PSICQUIC service to download data from IntAct ftp by selecting database argument "IntActFTP". This is much faster but larger requires larger download and is more computationally intensive for processing. As of 7.09.2017 "IntActFTP" provides access to DIP data, while "imex" doesn't. If database "IntActFTP" is chosen only MITAB2.7 is available and \code{format} is ignored
##' @inheritParams MItools fullInteractome
##' @param taxid1 character (1L), taxonomy id of the species which interaction participants should belong to, default is 9606 (which is human)
##' @param taxid2 character (1L), taxonomy id of the species which interaction participants should belong to, default is 9606 (which is human)
##' @return object of class `input class`_interSpeciesInteractome containing data.table containing molecular interaction data in either of these two formats:
##' @return if \code{clean} is TRUE: contains columns as described in \code{\link{cleanMITAB}};
##' @return if \code{clean} is FALSE: contains a standard set of columns for MITAB2.5 or MITAB2.7 depending on \code{format};
##' @seealso \code{\link{fullInteractome}}
##' @import data.table
##' @export interSpeciesInteractome
##' @export print.RAW_MItab25_interSpeciesInteractome
##' @export print.RAW_MItab27_interSpeciesInteractome
##' @export print.clean_MItab25_interSpeciesInteractome
##' @export print.clean_MItab27_interSpeciesInteractome
##' @examples
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.5 format, clean and select specific columns
##' interSpecies = interSpeciesInteractome(taxid1 = 9606,  taxid2 = 10239, database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE)
##'
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.7 format not using PSICQUIC (using IntAct ftp), clean and select specific columns; save it to the specific directory inside working directory
##' interSpecies = interSpeciesInteractome(taxid1 = 9606,  taxid2 = 10239, database = "IntActFTP", format = "tab27", clean = TRUE, protein_only = TRUE, directory = "./data/")
interSpeciesInteractome = function(MITABdata = NULL, taxid1 = 9606, taxid2 = 10239, database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE, directory = NULL, releaseORdate = NULL){
  # if the interaction data for species taxid and from database is not saved in the library - queryPSICQUIC for interaction data for taxid interactions in the database and in MITAB2.5 format, save results to the library
  if(is.null(MITABdata)){
    if(database == "IntActFTP"){
      if(is.null(directory)){
        pkg_dir = paste0(.libPaths(), "/MItools", "/data/")[1]
        # create data directory in /default.library/queryPSICQUIC/ if it doesn't exist
        if(!dir.exists(pkg_dir)) dir.create(pkg_dir)
        # find out last release date if the database is IntActFTP and releaseORdate = NULL, generate dir_last_release
        dir_last_release = generateDirName(database, releaseORdate, pkg_dir)
      } else {
        dir_last_release = generateDirName(database, releaseORdate, directory)
      }
      # create directory for the last release date
      if(is.null(releaseORdate)) {
        if(!dir.exists(dir_last_release)) dir.create(dir_last_release)
      } else {
        if(!dir.exists(dir_last_release)) stop(paste0("no data for IntAct release or date: ", releaseORdate," in the directory: ", directory))
      }
      full_interactome = loadIntActFTP(dir_last_release)
      taxids1 = loadTaxIDAllLower(taxid = taxid1, dir = dir_last_release)
      taxids1 = c(taxids$AllLower, taxids$input_taxid)
      taxids2 = loadTaxIDAllLower(taxid = taxid2, dir = dir_last_release)
      taxids2 = c(taxids$AllLower, taxids$input_taxid)
      full_interactome$data[, Taxid_interactor_A_clean := gsub("taxid:|\\(.*$","",`Taxid interactor A`)]
      full_interactome$data[, Taxid_interactor_B_clean := gsub("taxid:|\\(.*$","",`Taxid interactor B`)]
      full_interactome$data = all.intact$data[(Taxid_interactor_A_clean %in% taxids1 &
                                          Taxid_interactor_B_clean %in% taxids2) |
                                          (Taxid_interactor_A_clean %in% taxids2 &
                                             Taxid_interactor_B_clean %in% taxids1), ]
      full_interactome$data[, Taxid_interactor_A_clean := NULL]
      full_interactome$data[, Taxid_interactor_B_clean := NULL]
    } else {
    full_interactome = queryPSICQUICrlib(query = paste0("(taxidA:",taxid1," AND ", "taxidB:",taxid2, ") OR (taxidA:",taxid2," AND ", "taxidB:",taxid1, ")"),
                                         format = format,
                                         database = database,
                                         directory = directory,
                                         releaseORdate = releaseORdate)
    }
  }

  if(!is.null(MITABdata)) full_interactome = copy(MITABdata)

  # clean this data to make it more useble if clean is TRUE
  if(clean){
    full_interactome_clean = cleanMITAB(full_interactome)
    # filter out non-proteins
    if(protein_only){
      # removing interactions if at least one interactor has non-uniprot id
      full_interactome_clean$data = full_interactome_clean$data[interactor_IDs_databases_A == "uniprotkb" & interactor_IDs_databases_B == "uniprotkb",]
    }
    # reorder interactions so that all entities in IDs_interactor_A are of the same species and IDs_interactor_B are of the other
    full_interactome_clean$data
    full_interactome_clean$data = reorderMITAB27(full_interactome_clean$data)
    full_interactome_clean$taxid1 = taxid1
    full_interactome_clean$taxid2 = taxid2
    full_interactome_clean$protein_only = protein_only
    class(full_interactome_clean) = paste0(class(full_interactome_clean),"_interSpeciesInteractome")
    return(full_interactome_clean)
  }

  # if clean is FALSE return the interactome date in the raw MITAB format
  if(!clean){
    full_interactome$taxid = taxid
    full_interactome$protein_only = protein_only
    class(full_interactome) = paste0(class(full_interactome),"_interSpeciesInteractome")
    return(full_interactome)
  }

}

#print methods
print.RAW_MItab25_interSpeciesInteractome = function(data){
  cat(paste0("\n` Object of class RAW_MItab25_interSpeciesInteractome, contains interactions between molecular of taxid1: ", data$taxid1, " and taxid2:", data$taxid2," , \nnot ordered (the same species can be in IDs_interactor_A or IDs_interactor_B in different pairs), proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.RAW_MItab27_interSpeciesInteractome = function(data){
  cat(paste0("\n` Object of class RAW_MItab27_interSpeciesInteractome, contains interactions between molecular of taxid1: ", data$taxid1, " and taxid2:", data$taxid2," , \nnot ordered (the same species can be in IDs_interactor_A or IDs_interactor_B in different pairs), proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab25_interSpeciesInteractome = function(data){
  cat(paste0("\n` Object of class clean_MItab25_interSpeciesInteractome, contains interactions between molecular of taxid1: ", data$taxid1, " and taxid2:", data$taxid2," , proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_interSpeciesInteractome = function(data){
  cat(paste0("\n` Object of class clean_MItab27_interSpeciesInteractome, contains interactions between molecular of taxid1: ", data$taxid1, " and taxid2:", data$taxid2," , proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
