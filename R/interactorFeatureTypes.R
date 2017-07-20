##' discover types of interactor features in MITAB 2.7
##' @name interactorFeatureTypes
##' @author Vitalii Kleshchevnikov
##' @param mitab data.table containing molecular interaction data in MITAB 2.5 or 2.7 formats. Details: \code{\link{queryPSICQUIC}}
##' @details \code{interactorFeatureTypes} extracts types of interactor features from Features_interactor_A and Features_interactor_B in MITAB 2.7 format \link{https://psicquic.github.io/MITAB27Format.html}
##' @details All interactor features are decribed in the molecular interaction ontology (MI): \link{http://www.ebi.ac.uk/ols/ontologies/MI/terms?obo_id=MI:0116}
##' @import data.table
##' @importFrom gsubfn strapplyc
##' @export interactorFeatureTypes
##' @seealso \code{\link{cleanMITAB}}
interactorFeatureTypes = function(mitab){
  # discover types of features
  a = unique(unlist(mitab[, gsubfn::strapplyc(Features_interactor_A, "^([[:alnum:][:blank:]]*):|\\|([[:alnum:][:blank:]]*):", simplify = TRUE)]))
  b = unique(unlist(mitab[, gsubfn::strapplyc(Features_interactor_B, "^([[:alnum:][:blank:]]*):|\\|([[:alnum:][:blank:]]*):", simplify = TRUE)]))
  ab = unique(c(a,b))
  ab[order(ab)]
}
