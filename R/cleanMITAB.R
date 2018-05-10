##' Process molecular interaction data in MITAB 2.5 or 2.7 formats
##' @name cleanMITAB
##' @author Vitalii Kleshchevnikov
##' @description \code{cleanMITAB} extracts columns of molecular interaction data in MITAB 2.5 or 2.7 formats. For simple MITAB 2.5, columns are interactor Uniprot IDs(or other ID as specified by interactor_IDs_databases), interactor taxonomy IDs, Publication Identifiers, Confidence values and generates unique pair ID. MITAB 2.7 includes  additional information: see details.
##' @param MITABdata object of class "RAW_MItab25" or "RAW_MItab27" (list) containing molecular interaction data as returned by \code{\link{queryPSICQUICrlib}}
##' @return Object of S3 class "clean_MItab25" or "clean_MItab27" containing interaction data in a data.table. Within each interacting pair all interactor attributes are sorted according to alphanumerically sorted interactor UniprotKB accessions
##' @details Output column description (MITAB 2.5):
##' @details \code{pair_id} - unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs
##' @details \code{IDs_interactor_A}, \code{IDs_interactor_B} - interacting molecule ID
##' @details \code{interactor_IDs_databases_A}, \code{interactor_IDs_databases_B} - database that provides interacting molecule ID such as UniProt, ChEMBL or IntAct (IntAct: when the molecule cannot be mapped to ID in any other resource)
##' @details \code{Taxid_interactor_A}, \code{Taxid_interactor_B} - taxonomic species that interacting molecule belongs to
##' @details \code{Publication_Identifiers} - pubmed ID of the papers that has reported the interaction
##' @details \code{Confidence_values} - MIscore. Details: \link{https://psicquic.github.io/MITAB25Format.html}, \link{http://www.ebi.ac.uk/intact/pages/faq/faq.xhtml}, \link{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316181/}
##' @details if MITAB 2.7 then the following columns are added:
##' @details \code{Host_organisms} - taxid of organism in which the interaction detections experiment was performed
##' @details \code{bait_prey_status_A}, \code{bait_prey_status_B} - specifies the experimental role of a protein: a bait, a prey, a neutral component or else (\link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0495})
##' @details \code{Interaction_detection_methods} - Method to determine the interaction. \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0001}
##' @details \code{Interaction_types} - \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0190}
##' @details \code{Interaction_identifiers} - Interaction identifier given by database (e.g. intact:EBI-1547321|imex:IM-11865-3)
##' @details \code{Expansion_methods} - The method by which complex n-ary data is expanded into binary data. \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:1059}
##' @details \code{Features_interactor_A}, \code{Features_interactor_B} - Property of a subsequence that may interfere with the binding of a molecule. \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0116}
##' @details \code{Identification_method_participant_A}, \code{Identification_method_participant_B} - Method to determine the molecules involved in the interaction. http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0002
##' @details \code{binding_region_A_start}, \code{binding_region_A_end}, \code{binding_region_B_start}, \code{binding_region_B_end} - start and end position of binding-associated features
##' @details \code{binding_region_A_type}, \code{binding_region_B_type} - Type of the binding region. Details: \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0117}, \link{https://psicquic.github.io/MITAB27Format.html}, \code{\link{MITABregionFeature}}
##' @import data.table
##' @export cleanMITAB
##' @export print.clean_MItab25
##' @export print.clean_MItab27
##' @seealso \code{\link{interactorFeatureTypes}}
cleanMITAB = function(MITABdata){

  ###############################################################################
  # identifying MITAB format
  miformat = NA
  if(class(MITABdata) == "RAW_MItab25") {
    if(class(MITABdata$data)[1] == "data.table") if(ncol(MITABdata$data) == 15) miformat = "2.5"
  } else {
    if(class(MITABdata) == "RAW_MItab27") {
      if(class(MITABdata$data)[1] == "data.table") if(ncol(MITABdata$data) == 42) miformat = "2.7"
    } else {
      if(grepl("clean_MItab",class(cleanMITAB))) {
        message("the data has already been cleaned")
        return(MITABdata)
      } else stop("data is neither class RAW_MItab25 not RAW_MItab27")
    }
  }
  if(is.na(miformat)) stop("data is of the right class but not in MITAB 2.5 or 2.7 format, check if any columns were added or deleted from the original query output")

  MITABdata$data = copy(MITABdata$data)

  ###############################################################################
  # cleaning MITAB 2.5
  if(miformat == "2.5"){
    MITABdata$data = cleanMITAB25(MITABdata$data)
  }

  ###############################################################################
  # cleaning MITAB 2.7
  if(miformat == "2.7"){
    MITABdata$data = cleanMITAB27(MITABdata$data)
  }

  MITABdata$data = unique(MITABdata$data)
  if(miformat == "2.5") class(MITABdata) = "clean_MItab25"
  if(miformat == "2.7") class(MITABdata) = "clean_MItab27"
  return(MITABdata)
}

#print methods
print.clean_MItab25 = function(data){
  cat(paste0("\n` Object of class clean_MItab25 (molecular interaction data), for query, file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.clean_MItab27 = function(data){
  cat(paste0("\n` Object of class clean_MItab27 (molecular interaction data), for query, file, format, databases, date: `\n"))
  print(data$metadata)
  if("subsetByMethodDetails" %in% names(data)) printSubsetByMethodDetails(data)
  if("subsetByPMIDsDetails" %in% names(data)) printSubsetByPMIDsDetails(data)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
