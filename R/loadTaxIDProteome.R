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
  filename = paste0("UniProt_accessions_taxid",taxid,"_release",uniprot_release)
  if(!file.exists(filename)) download(paste0("http://www.uniprot.org/uniprot/?query=taxonomy:",taxid,"&format=tab&columns=id"), filename)
  proteins = unlist(fread(filename)[,1])
  names(proteins) = NULL
  list(proteins = proteins, taxid = taxid)
}

##' load all lower taxonomy ID
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
  filename = paste0("AllLowerTaxID_taxid",taxid,"_release",uniprot_release)
  if(!file.exists(filename)) download(paste0("http://www.uniprot.org/taxonomy/?query=ancestor:",taxid,"&format=tab"), filename)
  proteins = unlist(fread(filename)[,1])
  names(proteins) = NULL
  if(with_description){
    description = fread(filename)
    return(list(AllLower = proteins, input_taxid = taxid, description = description))
  } else {
    return(list(AllLower = proteins, input_taxid = taxid))
  }
}
