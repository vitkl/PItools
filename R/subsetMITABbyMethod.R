##' Filter molecular interactions by interaction detection methods and / or participant identification method
##' @name subsetMITABbyMethod
##' @author Vitalii Kleshchevnikov
##' @description filter molecular interaction data by \code{Interaction_detection_methods} and / or participant \code{Identification_method}
##' @param MITABdata object of any clean_MItab class (the class that contains "clean_MItab" in it's name and is initially produced by  \code{\link[MItools]{MITABdata}})
##' @param Interaction_detection_methods character or character vector, Molecular Interaction ontology code for Interaction detection method, such as "MI:0018", default option results in no filtering by Interaction detection method. Details \link{http://purl.obolibrary.org/obo/MI_0001}
##' @param Identification_method character or character vector, Molecular Interaction ontology code for Identification method, such as "MI:0433". Applies the same Identification method filter to both participant_A and participant_B. Default option results in no filtering by Identification method. Details \link{http://purl.obolibrary.org/obo/MI_0002}
##' @param Identification_method_participant_A character or character vector, Molecular Interaction ontology code for Identification method for participant A, such as "MI:0433", default option results in no filtering by Interaction detection method for this participant. Parameter is ignored if \code{Identification_method} specified. Details \link{http://purl.obolibrary.org/obo/MI_0002}
##' @param Identification_method_participant_B character or character vector, Molecular Interaction ontology code for Identification method for participant B, such as "MI:0433", default option results in no filtering by Interaction detection method for this participant. Parameter is ignored if \code{Identification_method} specified. Details \link{http://purl.obolibrary.org/obo/MI_0002}
##' @param ontology_directory where local copy of the MI ontology is stored (if file not found it will be downloaded)
##' @param ontology_date character, such as "2017_08_25", use specific version of ontology (local copy must be present), by default (NULL) download the latest version. The function doesn't require internet access if this argument is provided.
##' @param inverse_filter logical, inverse filtering criteria
##' @return object of the same class as \code{MITABdata}, a subset of \code{MITABdata} that contains interactions for \code{Interaction_detection_methods} and / or participant \code{Identification_method}
##' @import data.table
##' @seealso \code{\link[MItools]{subsetMITABbyPMIDs}}), \code{\link[MItools]{subsetMITABbyID}})
##' @export subsetMITABbyMethod
##' @examples
##' # Download all human interaction data
##' full = fullInteractome(taxid = "9606", database = "IntActFTP",
##'          clean = TRUE, protein_only = TRUE)
##'
##' # subset two-hybrid interactions
##' two_hybrid = subsetMITABbyMethod(MITABdata = full,
##'                Interaction_detection_methods = "MI:0018")
##'
##' # subset all interactions but two-hybrid
##' NOT_two_hybrid = subsetMITABbyMethod(MITABdata = full,
##'                    Interaction_detection_methods = "MI:0018", inverse_filter = T)
##'
##' # subset affinity purification - mass spectrometry interactions
##' AP_MS = subsetMITABbyMethod(MITABdata = full,
##'                Interaction_detection_methods = "MI:0004",  Identification_method = "MI:0433")
subsetMITABbyMethod = function(MITABdata, Interaction_detection_methods = NULL, Identification_method = NULL, Identification_method_participant_A = NULL, Identification_method_participant_B = NULL, ontology_directory = "./", ontology_date = NULL, inverse_filter = F) {
  if(!grepl("clean_MItab",class(MITABdata))) stop("MITABdata is not of class clean_MItab27 or related clean_MItab class")

  MIontology = loadMIontology(directory = ontology_directory, date = ontology_date)

  if(!is.null(Interaction_detection_methods)){
    all_Interaction_detection_methods = get_descendants(MIontology, roots = Interaction_detection_methods, exclude_roots = FALSE)
    ind = integer()
    for (method in all_Interaction_detection_methods) {
      ind = c(ind, grep(method,MITABdata$data$Interaction_detection_methods))
    }
    ind = unique(ind)
    if(inverse_filter){
      MITABdata$data = MITABdata$data[-ind,]
    } else MITABdata$data = MITABdata$data[ind,]
  }

  if(!is.null(Identification_method)){
    all_Identification_method = get_descendants(MIontology, roots = Identification_method, exclude_roots = FALSE)
    indA = integer()
    indB = integer()
    for (method in all_Identification_method) {
      indA = c(ind, grep(method,MITABdata$data$Identification_method_participant_A))
      indB = c(ind, grep(method,MITABdata$data$Identification_method_participant_B))
    }
    ind2 = intersect(indA, indB)
    if(inverse_filter){
      MITABdata$data = MITABdata$data[-ind2,]
    } else MITABdata$data = MITABdata$data[ind2,]
  } else {
    if(!is.null(Identification_method_participant_A)){
      all_Identification_method_participant_A = get_descendants(MIontology, roots = Identification_method_participant_A, exclude_roots = FALSE)
      indA = integer()
      for (method in all_Identification_method_participant_A) {
        indA = c(ind, grep(method, MITABdata$data$Identification_method_participant_A))
      }
      if(inverse_filter){
        MITABdata$data = MITABdata$data[-indA,]
      } else MITABdata$data = MITABdata$data[indA,]
    }
    if(!is.null(Identification_method_participant_B)){
      all_Identification_method_participant_B = get_descendants(MIontology, roots = Identification_method_participant_B, exclude_roots = FALSE)
      indB = integer()
      for (method in all_Identification_method_participant_B) {
        indB = c(ind, grep(method, MITABdata$data$Identification_method_participant_B))
      }
      if(inverse_filter){
        MITABdata$data = MITABdata$data[-indB,]
      } else MITABdata$data = MITABdata$data[indB,]
    }
  }

  MITABdata$subsetByMethodDetails = list(Interaction_detection_methods = Interaction_detection_methods,
                                         Identification_method = Identification_method,
                                         Identification_method_participant_A = Identification_method_participant_A,
                                         Identification_method_participant_B = Identification_method_participant_B,
                                         ontology_directory = ontology_directory,
                                         ontology_date = ontology_date,
                                         inverse_filter = inverse_filter)
  return(MITABdata)

}

printSubsetByMethodDetails = function(data){
  cat(paste0("\n` Data was also filtered by Interaction_detection_methods and / or participant Identification_methods `\n"))
  print(data$subsetByMethodDetails[!sapply(data$subsetByMethodDetails, is.null)])
}
