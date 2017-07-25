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
    mitab = MITABregionFeature(mitab)

    # reorder by all interactor attribute columns by pair_id (alphanumeric order)
    mitab[, c("IDs_A_order", "IDs_B_order") := tstrsplit(pair_id, "\\|")]
    mitab[IDs_interactor_A == IDs_B_order & IDs_interactor_B == IDs_A_order,
          c("IDs_interactor_A", "IDs_interactor_B",
            "interactor_IDs_databases_A", "interactor_IDs_databases_B",
            "Taxid_interactor_A", "Taxid_interactor_B",
            "bait_prey_status_A", "bait_prey_status_B",
            "Features_interactor_A", "Features_interactor_B",
            "Identification_method_participant_A", "Identification_method_participant_B",
            "binding_region_A_start", "binding_region_A_end", "binding_region_B_start", "binding_region_B_end",
            "binding_region_A_type", "binding_region_B_type") :=
            .(IDs_interactor_B, IDs_interactor_A,
               interactor_IDs_databases_B, interactor_IDs_databases_A,
               Taxid_interactor_B, Taxid_interactor_A,
               bait_prey_status_B, bait_prey_status_A,
               Features_interactor_B, Features_interactor_A,
               Identification_method_participant_B, Identification_method_participant_A,
               binding_region_B_start, binding_region_B_end, binding_region_A_start, binding_region_A_end,
               binding_region_B_type, binding_region_A_type)]

    # start by keeping only relevant columns
    mitab = unique(mitab[, .(IDs_interactor_A, IDs_interactor_B,
                             interactor_IDs_databases_A, interactor_IDs_databases_B,
                             Taxid_interactor_A, Taxid_interactor_B,
                             Publication_Identifiers, Confidence_values,
                             Host_organisms,
                             bait_prey_status_A, bait_prey_status_B,
                             Interaction_detection_methods, Interaction_types, Interaction_identifiers, Expansion_methods,
                             Features_interactor_A, Features_interactor_B,
                             Identification_method_participant_A, Identification_method_participant_B,
                             binding_region_A_start, binding_region_A_end, binding_region_B_start, binding_region_B_end,
                             binding_region_A_type, binding_region_B_type,
                             pair_id)])
    mitab = MITABregionFeature(mitab)
  }
}
