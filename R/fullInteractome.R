##' Retrieve interactome (proteins only or not) of a given taxonomic species
##' @name fullInteractome
##' @author Vitalii Kleshchevnikov
##' @description Retrieve interactome (proteins only or not) of a given taxonomic species from a specified database. Interactome can be additionally cleaned and includes only specific information: \code{\link{cleanMITAB}}
##' @details \code{taxid} is used to query specified database using PSICQUIC client, only interactions in which both participants belong the \code{taxid} are retured (\code{"taxidA:9606 AND taxidB:9606"}, not \code{"species:9606"}). Details: \code{\link{queryPSICQUIC}}
##' @details \code{fullInteractome} can be used to retrive interactome data using PSICQUIC service using \code{\link{queryPSICQUIC}}, clean and select specific columns using \code{\link{cleanMITAB}} and filter resulting dataset for protein-protein interaction only. This is the default option.
##' @details Alternatively, \code{fullInteractome} can only retrive interactome data using PSICQUIC service without cleaning of filtering.
##' @details Another option is to supply \code{MITABdata} to be cleaned and filtered
##' @details Finally, you can avoid using PSICQUIC service and download data from IntAct ftp by selecting database argument "IntActFTP". This is much faster but larger requires larger download (>3Gb) and is more computationally intensive for processing (clean = T). As of 7.09.2017 "IntActFTP" provides access to DIP data, while "imex" doesn't.
##' @param MITABdata object of class "RAW_MItab25" or "RAW_MItab27" (list) containing molecular interaction data as returned by \code{\link{queryPSICQUICrlib}} or \code{\link{loadIntActFTP}}, default in NULL
##' @param taxid character (1L), taxonomy id of the species which interaction participants should belong to, default is "9606" (which is human)
##' @param database character (1L), argument for \code{\link{queryPSICQUIC}}, PSICQUIC-compliant database to query for interactions. The default is "imex" alternative to which is "IntActFTP"
##' @param format character (1L), argument for \code{\link{queryPSICQUIC}}, default is "tab25"
##' @param clean logical (1L), if TRUE extract specific information using \code{\link{cleanMITAB}}, default is TRUE
##' @param protein_only logical (1L), if TRUE the interaction participants are restricted to proteins (exclude other types of molecules such as RNA or small molecules), default is TRUE
##' @param directory directory where to store the data, if NULL the data is stored in <R-package-library>/PItools/data
##' @param releaseORdate character, if data has already been downloaded: which IntAct release or download date to read
##' @param remove_obsolete_id logical (1L), remove interactions in which one of the partners is encoded as obsolete UniProtKB accession (ID), not implemented properly: will never finish.
##' @param within_species logical (1L), return interactions only between proteins of \code{taxid}. If FALSE returns interactions for \code{taxid} proteins both within \code{taxid} and with proteins from other species
##' @return object of class `input class`_fullInteractome containing data.table containing molecular interaction data in either of these two formats:
##' @return if \code{clean} is TRUE: contains columns as described in \code{\link{cleanMITAB}};
##' @return if \code{clean} is FALSE: contains a standard set of columns for MITAB2.5 or MITAB2.7 depending on \code{format};
##' @seealso \code{\link{interSpeciesInteractome}}
##' @import data.table
##' @export fullInteractome
##' @export print.RAW_MItab25_fullInteractome
##' @export print.RAW_MItab27_fullInteractome
##' @export print.clean_MItab25_fullInteractome
##' @export print.clean_MItab27_fullInteractome
##' @examples
##' {
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.5 format, clean and select specific columns (commented because takes a lot of time compared to IntActFTP)
##' # full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE)
##'
##' # Do the same using IntAct ftp instead of the PSICQUIC webservice, except that, IntAct ftp option always outputs tab27 format
##' full = fullInteractome(taxid = "9606", database = "IntActFTP", clean = TRUE, protein_only = TRUE)
##'
##' # retrive a full set of human (9606) protein-protein interactions from IMEx databases in MITAB2.5 format, clean and select specific columns; save it to the specific directory inside working directory
##' # full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE, directory = "./data/")
##' }
fullInteractome = function(MITABdata = NULL, taxid = 9606, database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE, directory = NULL, releaseORdate = NULL, remove_obsolete_id = F, within_species = T){
  # if the interaction data for species taxid and from database is not saved in the library - queryPSICQUIC for interaction data for taxid interactions in the database and in MITAB2.5 format, save results to the library
  if(database == "IntActFTP"){
    if(is.null(directory)){
      pkg_dir = paste0(.libPaths(), "/PItools", "/data/")[1]
      # create /data/ directory in /default.library/queryPSICQUIC/ if it doesn't exist
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
      if(!dir.exists(dir_last_release)) stop(paste0("no data for IntAct release or date: ", releaseORdate," in the directory: ", directory, ", set release to NULL do load the latest release"))
    }
    if(!is.null(MITABdata)) full_interactome = copy(MITABdata) else full_interactome = loadIntActFTP(dir_last_release)

    taxids = loadTaxIDAllLower(taxid = taxid, dir = dir_last_release)
    taxids = c(taxids$AllLower, taxids$input_taxid)
    full_interactome$data[, Taxid_interactor_A_clean := gsub("taxid:|\\(.*$","",`Taxid interactor A`)]
    full_interactome$data[, Taxid_interactor_B_clean := gsub("taxid:|\\(.*$","",`Taxid interactor B`)]
    if(within_species){
      query = paste0("taxidA:",taxid," AND ", "taxidB:",taxid)
      full_interactome$data = full_interactome$data[Taxid_interactor_A_clean %in% taxids & Taxid_interactor_B_clean %in% taxids, ]
    } else {
      query = paste0("taxidA:",taxid," OR ", "taxidB:",taxid)
      full_interactome$data = full_interactome$data[Taxid_interactor_A_clean %in% taxids | Taxid_interactor_B_clean %in% taxids, ]
    }

    full_interactome$data[, Taxid_interactor_A_clean := NULL]
    full_interactome$data[, Taxid_interactor_B_clean := NULL]
    full_interactome$metadata = data.table(query = query,
                                           file = paste0(dir_last_release,"intact",gsub("^.*IntActRelease_|/","", dir_last_release),".txt.gz"),
                                           format = "tab27",
                                           all.databases = paste0("IntActFTP_", names(table(full_interactome$data$`Source database(s)`))),
                                           n.interactions.in.database = table(full_interactome$data$`Source database(s)`),
                                           database.not.active = "NULL")
  } else {
    if(!is.null(MITABdata)) {
      full_interactome = copy(MITABdata)
      taxids = loadTaxIDAllLower(taxid = taxid, dir = dir_last_release)
      taxids = c(taxids$AllLower, taxids$input_taxid)
      full_interactome$data[, Taxid_interactor_A_clean := gsub("taxid:|\\(.*$","",`Taxid interactor A`)]
      full_interactome$data[, Taxid_interactor_B_clean := gsub("taxid:|\\(.*$","",`Taxid interactor B`)]
      if(within_species){
        query = paste0("taxidA:",taxid," AND ", "taxidB:",taxid)
        full_interactome$data = full_interactome$data[Taxid_interactor_A_clean %in% taxids & Taxid_interactor_B_clean %in% taxids, ]
      } else {
        query = paste0("taxidA:",taxid," OR ", "taxidB:",taxid)
        full_interactome$data = full_interactome$data[Taxid_interactor_A_clean %in% taxids | Taxid_interactor_B_clean %in% taxids, ]
      }

      full_interactome$data[, Taxid_interactor_A_clean := NULL]
      full_interactome$data[, Taxid_interactor_B_clean := NULL]
      full_interactome$metadata$query = query
    } else {
      if(within_species){
        query = paste0("taxidA:",taxid," AND ", "taxidB:",taxid)
      } else {
        query = paste0("taxidA:",taxid," OR ", "taxidB:",taxid)
      }
      full_interactome = queryPSICQUICrlib(query = query,
                                           format = format,
                                           database = database,
                                           directory = directory,
                                           releaseORdate = releaseORdate)
    }
  }

  # clean this data to make it more useble if clean is TRUE
  if(clean){
    full_interactome_clean = cleanMITAB(full_interactome)
    # filter out non-proteins
    if(protein_only){
      # removing interactions if at least one interactor has non-uniprot id
      full_interactome_clean$data = full_interactome_clean$data[interactor_IDs_databases_A == "uniprotkb" & interactor_IDs_databases_B == "uniprotkb",]
    }
    full_interactome_clean$taxid = taxid
    full_interactome_clean$protein_only = protein_only
    full_interactome_clean$remove_obsolete_id = remove_obsolete_id
    class(full_interactome_clean) = paste0(class(full_interactome_clean),"_fullInteractome")
    if(remove_obsolete_id){
      full_interactome_clean = removeInteractionObsoleteID(full_interactome_clean, dir = directory)
    }
    return(full_interactome_clean)
  } else {   # if clean is FALSE return the interactome date in the raw MITAB format
    full_interactome$taxid = taxid
    full_interactome$protein_only = protein_only
    class(full_interactome) = paste0(class(full_interactome),"_fullInteractome")
    return(full_interactome)
  }

}

#print methods
print.RAW_MItab25_fullInteractome = function(data){
  cat(paste0("\n` Object of class RAW_MItab25_fullInteractome, for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.RAW_MItab27_fullInteractome = function(data){
  cat(paste0("\n` Object of class RAW_MItab27_fullInteractome, for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab25_fullInteractome = function(data){
  cat(paste0("\n` Object of class clean_MItab25_fullInteractome, for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27_fullInteractome = function(data){
  cat(paste0("\n` Object of class clean_MItab27_fullInteractome, for taxid: ", data$taxid, ", proteins only: ", data$protein_only," `\n"))
  cat(paste0("\n` file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
