##' load the list of all proteins for given taxonomy ID (any level)
##' @name loadTaxIDProteins
##' @author Vitalii Kleshchevnikov
##' @description loadTaxIDProteins function loads all UniProt accessions for a given taxonomy ID
##' @param dir directory where to save/look for the local copy
##' @return \code{loadTaxIDProteins} returns list containing the character vector of all UniProt accessions for a given taxonomy ID (saves it to a file) and the taxonomy ID (taxid)
##' @import data.table
##' @importFrom httr GET
##' @importFrom httr headers
##' @importFrom downloader download
##' @export loadTaxIDProteins
##' @examples loadTaxIDProteins(taxid = 9606, dir = "./")
##' @author Vitalii Kleshchevnikov
loadTaxIDProteins = function(taxid, dir){
  uniprot_release = httr::headers(httr::GET("http://www.uniprot.org/"))$`x-uniprot-release`
  filename = paste0(dir, "UniProt_accessions_taxid",taxid,"_release",uniprot_release)
  if(!file.exists(filename)) download(paste0("http://www.uniprot.org/uniprot/?query=taxonomy:",taxid,"&format=tab&columns=id"), filename)
  proteins = unlist(fread(filename)[,1])
  names(proteins) = NULL
  list(proteins = proteins, taxid = taxid)
}

##' load all lower taxonomy ID (all descendants)
##' @name loadTaxIDAllLower
##' @author Vitalii Kleshchevnikov
##' @description \code{loadTaxIDAllLower} loads all load all lower taxonomy ID for a given taxonomy ID
##' @param file directory where to save/look for the local copy
##' @return \code{loadTaxIDAllLower} returns list containing the character vector of all lower taxonomy ID for a given taxonomy ID (saves it to a file) and the input taxonomy ID (taxid) and optionally, a data.table containing the description of all lower
##' @import data.table
##' @importFrom httr GET
##' @importFrom httr headers
##' @importFrom downloader download
##' @export loadTaxIDAllLower
##' @examples loadTaxIDAllLower(taxid = 9606, dir = "./")
##' loadTaxIDAllLower(taxid = 9606, dir = "./", with_description = T)
##' @author Vitalii Kleshchevnikov
loadTaxIDAllLower = function(taxid, dir, with_description = F){
  uniprot_release = httr::headers(httr::GET("http://www.uniprot.org/"))$`x-uniprot-release`
  filename = paste0(dir,"AllLowerTaxID_taxid",taxid,"_release",uniprot_release)
  if(!file.exists(filename)) download(paste0("http://www.uniprot.org/taxonomy/?query=ancestor:",taxid,"&format=tab"), filename)
  if(file.size(filename) == 0) taxids = NULL else {
    taxids = unlist(fread(filename)[,1])
    names(taxids) = NULL
  }
  if(with_description){
    description = fread(filename)
    return(list(AllLower = taxids, input_taxid = taxid, description = description))
  } else {
    return(list(AllLower = taxids, input_taxid = taxid))
  }
}

##' load the list of all protein ID in UniProtKB
##' @name loadAllIDProteins
##' @author Vitalii Kleshchevnikov
##' @description loadAllIDProteins function loads all UniProt accessions for a given taxonomy ID
##' @param dir directory where to save/look for the local copy
##' @return \code{loadAllIDProteins} returns list containing the character vector of all UniProt accessions (saves it to a file)
##' @import data.table
##' @importFrom httr GET
##' @importFrom httr headers
##' @importFrom downloader download
##' @export loadAllIDProteins
##' @examples loadAllIDProteins(dir = "./")
##' @author Vitalii Kleshchevnikov
loadAllIDProteins = function(dir){
  uniprot_release = httr::headers(httr::GET("http://www.uniprot.org/uniprot/?query=*"))$`x-uniprot-release`
  total_results = as.numeric(httr::headers(httr::GET("http://www.uniprot.org/uniprot/?query=*"))$`x-total-results`)
  filename = paste0(dir, "All_UniProt_accessions_release",uniprot_release, ".txt")
  if(!file.exists(filename) | file.size(filename) == 0){
    system(paste0("touch ", filename))
    for (i in seq(1, total_results, 1e6)) {
      query = paste0("http://www.uniprot.org/uniprot/?format=list&limit=1000000&offset=",i)
      temp_file = paste0(tempdir(), "All_UniProt_accessions",i, ".txt")
      download.file(query, temp_file)
      system(paste0("cat ",temp_file," >> ",filename))
    }
  }
  #proteins = fread(filename, stringsAsFactors = F, header = F)
  filename
}

##' Remove interactions if one participant has obsolete Uniprot accession (ID)
##' @name removeInteractionObsoleteID
##' @author Vitalii Kleshchevnikov
##' @description filter molecular interaction data by Pubmed ID of publications
##' @param MITABdata object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by \code{\link[MItools]{cleanMITAB}})
##' @param dir directory where to save/look for the local copy
##' @import data.table
##' @export removeInteractionObsoleteID
removeInteractionObsoleteID = function(MITABdata, dir = "./"){
  interactors = extractInteractors(MITABdata)
  # shave off isoform ids
  interactors = gsub("-[[:digit:]]+$", "", interactors)
  # shave off post-processed chain ids
  interactors = gsub("-PRO_[[:digit:]]+$", "", interactors)

  filename = loadAllIDProteins(dir = dir)

  which_interactors = sapply(interactors, function(interactor, filename){
    length(system(paste("grep -w", interactor, filename), intern = T)) >=1
  }, filename)
  interactors = interactors[which_interactors]
  MITABdata = subsetMITABbyID(MITABdata, interactors, within_seed = T)
  MITABdata
}

##' Remove interactions if one participant doesn't have FASTA sequence
##' @name removeInteractionNoFASTA
##' @author Vitalii Kleshchevnikov
##' @description filter molecular interaction data by the list of FASTA sequences in AAStringSet
##' @param MITABdata object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by \code{\link[MItools]{cleanMITAB}})
##' @param fasta AAStringSet, names are UniProtKB accessions
##' @import data.table
##' @import Biostrings
##' @export removeInteractionNoFASTA
removeInteractionNoFASTA = function(MITABdata, fasta){
  interactors = extractInteractors(MITABdata)
  has_interactors = interactors %in% names(fasta)
  interactors = interactors[has_interactors]
  MITABdata = subsetMITABbyID(MITABdata, interactors, within_seed = T)
  MITABdata
}
