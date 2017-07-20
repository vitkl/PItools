##' \code{cleanMITAB} extracts and cleans a set of columns of MITAB2.7
##' @description function called internally by \code{\link{cleanMITAB}} if the format is MITAB2.7
##' @param mitab data.table containing molecular interaction data in MITAB 2.7 format
##' @author Vitalii Kleshchevnikov
##' @import data.table
cleanMITAB27 = function(mitab){
  # changing column names to data.table-compatible format
  {
    colnames(mitab) = gsub(" ","_",colnames(mitab))
    colnames(mitab) = gsub("\\(|\\)","",colnames(mitab))
    colnames(mitab) = gsub("#","",colnames(mitab))
  }
  # cleaning Taxid "taxid:9606(human)|taxid:9606(Homo sapiens)" to 9606
  {
    mitab[, Taxid_interactor_A := gsub("taxid:|\\(.*$","",Taxid_interactor_A)]
    mitab[, Taxid_interactor_B := gsub("taxid:|\\(.*$","",Taxid_interactor_B)]
    mitab[, Host_organisms := gsub("taxid:|\\(.*$","",Host_organisms)]
    # saving identifier types and cleaning interactor ids
    mitab[, interactor_IDs_databases_A := gsub(":.*$","",IDs_interactor_A)]
    mitab[, interactor_IDs_databases_B := gsub(":.*$","",IDs_interactor_B)]
    mitab[, IDs_interactor_A := gsub("^.*:","",IDs_interactor_A)]
    mitab[, IDs_interactor_B := gsub("^.*:","",IDs_interactor_B)]
    # isoform "-1" is a canonical sequence, IntAct uses isoform "-1" when it's clear that the isoform is "-1" and a canonical identifier if it's not clear which isoform was used in the experiment. Removing isoform sign "-1":
    mitab[, IDs_interactor_A := gsub("-1$", "", IDs_interactor_A)]
    mitab[, IDs_interactor_B := gsub("-1$", "", IDs_interactor_B)]
    # cleaning other information
    mitab[, bait_prey_status_A := gsub("^.*\\(|\\)","",Experimental_roles_interactor_A)]
    mitab[, bait_prey_status_B := gsub("^.*\\(|\\)","",Experimental_roles_interactor_B)]
    mitab[, Publication_Identifiers := gsub("^.*pubmed:|\\|.*$","",Publication_Identifiers)]
    mitab[, Confidence_values := gsub("^intact-miscore:","",Confidence_values)]
    mitab[, Confidence_values := gsub("-","NA",Confidence_values)]
    # supress expected warning (NA introduced by coersion to numeric) to avoid confusion
    suppressWarnings({mitab[, Confidence_values := as.numeric(Confidence_values)]})
    #mitab[, Interaction_identifiers := unlist(gsubfn::strapplyc(Interaction_identifiers,"EBI-[[:digit:]]+",simplify = T)), by =Interaction_identifiers]
    # generating unique identifier for interacting pairs
    mitab[, pair_id := apply(data.table(IDs_interactor_A,IDs_interactor_B,stringsAsFactors = F), 1,
                             function(a) { z = sort(a)
                             paste0(z[1],"|",z[2]) })]
  }
  {
    # extract region sufficient to interact or mutation affecting interaction information from Features_interactor_A and Features_interactor_B
    # start by keeping only relevant columns
    mitab = unique(mitab[, .(pair_id,
                             IDs_interactor_A, IDs_interactor_B,
                             interactor_IDs_databases_A, interactor_IDs_databases_B,
                             Taxid_interactor_A, Taxid_interactor_B,
                             Publication_Identifiers, Confidence_values,
                             Host_organisms,
                             bait_prey_status_A, bait_prey_status_B,
                             Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                             Features_interactor_A, Features_interactor_B,
                             Identification_method_participant_A, Identification_method_participant_B)])
    mitab = MITABregionFeature(mitab)
  }
}
