##' load Molecular Interaction ontology
##' @name loadMIontology
##' @author Vitalii Kleshchevnikov
##' @description download and read Molecular Interaction ontology into R using \code{ontologyIndex} package
##' @param directory where to keep local copy of the MI ontology
##' @param date character, such as "2017_08_25", use specific version of ontology (local copy must be present), by default (NULL) download the latest version
##' @param propagate_relationships parameter for reading ontology into R, details: \code{\link[ontologyIndex]{get_ontology}}
##' @param extract_tags parameter for reading ontology into R, details: \code{\link[ontologyIndex]{get_ontology}}
##' @return \code{MITABdata}
##' @import ontologyIndex
##' @importFrom jsonlite fromJSON
##' @importFrom downloader download
##' @export loadMIontology
loadMIontology = function(directory = "./", date = NULL, propagate_relationships = "is_a", extract_tags = "minimal") {
  if(is.null(date)){
    version = jsonlite::fromJSON("https://www.ebi.ac.uk/ols/api/ontologies/mi")$config$annotations$date
    version = as.POSIXct(strptime(version, "%d:%m:%Y %H:%M"))
    filename = paste0(directory, gsub("-", "_", as.Date(version)),"_mi.obo")
  } else {
    filename = paste0(directory, date,"_mi.obo")
  }

  if(!file.exists(filename)) {
    message("downloading MI ontology from http://purl.obolibrary.org/obo/mi.obo")
    downloader::download("http://purl.obolibrary.org/obo/mi.obo", filename)
  } else message("loading local copy of MI ontology")
  MIontology = get_ontology(filename, propagate_relationships = propagate_relationships,
                            extract_tags = extract_tags)
  return(MIontology)
}

##' browse Molecular Interaction ontology
##' @name browseMIontology
##' @author Vitalii Kleshchevnikov
##' @export browseMIontology
browseMIontology = function(){
  browseURL(url = "https://www.ebi.ac.uk/ols/ontologies/mi")
}
