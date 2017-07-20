##' \code{cleanMITAB} extracts and cleans a set of columns of MITAB2.5
##' @description function called internally by \code{\link{cleanMITAB}} if the format is MITAB2.5. No column names.
##' @param mitab data.table containing molecular interaction data in MITAB 2.5 format
##' @author Vitalii Kleshchevnikov
##' @import data.table
cleanMITAB25 = function(mitab){
  # cleaning Taxid "taxid:9606(human)|taxid:9606(Homo sapiens)" to 9606
  mitab[, Taxid_interactor_A := gsub("taxid:|\\(.*$","",V10)]
  mitab[, Taxid_interactor_B := gsub("taxid:|\\(.*$","",V11)]
  # saving identifier types and cleaning interactor ids
  mitab[, interactor_IDs_databases_A := gsub(":.*$","",V1)]
  mitab[, interactor_IDs_databases_B := gsub(":.*$","",V2)]
  mitab[, IDs_interactor_A := gsub("^.*:","",V1)]
  mitab[, IDs_interactor_B := gsub("^.*:","",V2)]
  # isoform "-1" is a canonical sequence, IntAct uses isoform "-1" when it's clear that the isoform is "-1" and a canonical identifier if it's not clear which isoform was used in the experiment. Removing isoform sign "-1":
  mitab[, IDs_interactor_A := gsub("-1$", "", IDs_interactor_A)]
  mitab[, IDs_interactor_B := gsub("-1$", "", IDs_interactor_B)]
  # cleaning other information
  mitab[, Publication_Identifiers := gsub("^.*pubmed:|\\|.*$","",V9)]
  mitab[, Confidence_values := gsub(".*intact-miscore:","",V15)]
  mitab[, Confidence_values := gsub("-","NA",Confidence_values)]
  # supress expected warning (NA introduced by coersion to numeric) to avoid confusion
  suppressWarnings({mitab[, Confidence_values := as.numeric(Confidence_values)]})
  #mitab[, Interaction_identifiers := unlist(gsubfn::strapplyc(Interaction_identifiers,"EBI-[[:digit:]]+",simplify = T)), by =Interaction_identifiers]
  # generating unique identifier for interacting pairs
  mitab[, pair_id := apply(data.table(IDs_interactor_A,IDs_interactor_B,stringsAsFactors = F), 1,
                           function(a) { z = sort(a)
                           paste0(z[1],"|",z[2]) })]
  mitab = mitab[,.(pair_id,
                   IDs_interactor_A, IDs_interactor_B,
                   interactor_IDs_databases_A, interactor_IDs_databases_B,
                   Taxid_interactor_A, Taxid_interactor_B,
                   Publication_Identifiers, Confidence_values)]
}
