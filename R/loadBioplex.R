##' load Bioplex data
##' @name loadBioplex
##' @author Vitalii Kleshchevnikov
##' @description \code{loadBioplex} saves and/or loads the data stored on Bioplex website \link{http://bioplex.hms.harvard.edu/downloadInteractions.php}. Depends on the tsv file containing the following columns:
##' Gene A	Gene B	Uniprot A	Uniprot B	Symbol A	Symbol B	p(Wrong)	p(No Interaction)	p(Interaction)
##' @param dir directory where to save/look for the local copy
##' @param uniprot_id If TRUE uses UniprotKB accessions to identify interacting partners. If FALSE uses gene id (ENTREZ gene id) to identify interacting partners.
##' @return \code{loadBioplex} saves loadBioplex data to a file, returns the object of class clean_MItab25
##' @import data.table
##' @importFrom R.utils gzip
##' @importFrom R.utils gunzip
##' @importFrom downloader download
##' @export loadBioplex
##' @examples loadBioplex("./")
##' @author Vitalii Kleshchevnikov
loadBioplex = function(dir = "./", url = "http://bioplex.hms.harvard.edu/data/BioPlex_2.3_interactionList.tsv", uniprot_id = T) {

  file_name = paste0(dir,gsub(".+/","","http://bioplex.hms.harvard.edu/data/BioPlex_2.3_interactionList.tsv"))
  file_name.gz = paste0(file_name, ".gz")
  # download Bioplex data
  if(!file.exists(file_name.gz)) {
    message("... dowloading from http://bioplex.hms.harvard.edu/ ...")
    downloader::download(url = url, destfile = file_name)
    R.utils::gzip(filename = file_name, destname = file_name.gz, remove = T, overwrite = T)
  } else {message("... loading local copy ...")}
  R.utils::gunzip(filename = file_name.gz, destname = file_name, remove = F, overwrite = T)
  bioplex = fread(file_name, header = T, stringsAsFactors = F)
  unlink(file_name)

  if(uniprot_id){
    bioplex = bioplex[,.(IDs_interactor_A = UniprotA, IDs_interactor_B = UniprotB,
                         interactor_IDs_databases_A = "uniprotkb", interactor_IDs_databases_B = "uniprotkb",
                         Taxid_interactor_A = 9606, Taxid_interactor_B = 9606,
                         Publication_Identifiers = paste("Bioplex", "GeneA", GeneA, "GeneB", GeneB, "SymbolA", SymbolA, "SymbolB", SymbolB, "p(Wrong)", `p(Wrong)`, "p(No Interaction)", `p(No Interaction)`, sep = "|"),
                         Confidence_values = `p(Interaction)`)]
  } else {
    bioplex = bioplex[,.(IDs_interactor_A = GeneA, IDs_interactor_B = GeneB,
                         interactor_IDs_databases_A = "entrez", interactor_IDs_databases_B = "entrez",
                         Taxid_interactor_A = 9606, Taxid_interactor_B = 9606,
                         Publication_Identifiers = paste("Bioplex", "UniprotA", UniprotA, "UniprotB", UniprotB, "SymbolA", SymbolA, "SymbolB", SymbolB, "p(Wrong)", `p(Wrong)`, "p(No Interaction)", `p(No Interaction)`, sep = "|"),
                         Confidence_values = `p(Interaction)`)]
  }


  # adding empty columns to fit the MItab27 format
  bioplex[, c("Host_organisms", "bait_prey_status_A", "bait_prey_status_B", "Interaction_detection_methods", "Interaction_types", "Interaction_identifiers", "Expansion_methods", "Features_interactor_A", "Features_interactor_B", "Identification_method_participant_A", "Identification_method_participant_B", "binding_region_A_start", "binding_region_A_end", "binding_region_B_start", "binding_region_B_end", "binding_region_A_type", "binding_region_B_type") := .(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)]

  # generating unique identifier for interacting pairs
  bioplex[, pair_id := apply(data.table(IDs_interactor_A,IDs_interactor_B,stringsAsFactors = F), 1,
                           function(a) { z = sort(a)
                           paste0(z[1],"|",z[2]) })]

  # reorder by all interactor attribute columns by pair_id (alphanumeric order)
  bioplex[, c("IDs_A_order", "IDs_B_order") := tstrsplit(pair_id, "\\|")]
  bioplex = reorderMITAB27(bioplex)

  MITABdata = list(data = bioplex, metadata = paste0("Summary of protein interactions in the BioPlex network following reanalysis of all data from BioPlex 2.0 plus an additional ~1500 unpublished AP-MS experiments: ", url))
  class(MITABdata) = "clean_MItab27"
  MITABdata
}
