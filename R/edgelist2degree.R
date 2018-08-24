##' find the degree for each protein given cleaned molecular interaction data
##' @name edgelist2degree
##' @author Vitalii Kleshchevnikov
##' @details \code{edgelist2degree} uses "IDs_interactor_A", "IDs_interactor_B" columns of a supplied data.table \code{mitab} to compute the degree. \code{mitab} can be generated using \code{\link{cleanMITAB}}
##' @details \code{edgelist2degree_slow} uses pair_id to compute degree. Unique pair_id is the unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs \code{\link{cleanMITAB}}
##' @param mitab data.table containing molecular interaction data, including "IDs_interactor_A" and "IDs_interactor_B" or "pair_id"
##' @param sep how are individual interactors in the pair_id reparated? default: "\\|"
##' @return sorted data.table containing degree for each of the nodes in the input, details: \code{mitab}
##' @import data.table
##' @export edgelist2degree
##' @export edgelist2degree_slow
##' @usage
##' edgelist2degree(mitab)
##' edgelist2degree_slow(mitab, sep = "\\|")
##' @examples
##' edgelist = data.table(pair_id = paste0(rep(c("P", "Q", "R"), each = 30),
##'                       "|",
##'                       sample(rep(c("Z", "X", "F"), each = 30)), 1:90))
##' edgelist2degree_slow(edgelist)
##'
##' # download full human interactome (clean = TRUE is necessary to produce the right input for edgelist2degree)
##' full = fullInteractome(taxid = "9606", database = "IntActFTP", format = "tab25", clean = TRUE, protein_only = TRUE)
##' degree = edgelist2degree(full$data)
edgelist2degree_slow = function(mitab, sep = "\\|", order_by_degree = T){
  if("pair_id" %in% colnames(mitab)){
    mitab_t = unique(copy(mitab[,.(pair_id)]))
    mitab_t[, c("ida", "idb") := tstrsplit(pair_id, sep)]
    degrees = mitab_t[, .(c(ida, idb))][, .N, by = V1]
    setnames(degrees, colnames(degrees), c("ID", "N"))
  } else stop("pair_id column not supplied")
  if(order_by_degree) setorder(degrees, N, ID)
  return(degrees)
}

edgelist2degree = function(mitab, order_by_degree = T){
  if(mean(c("IDs_interactor_A", "IDs_interactor_B") %in% colnames(mitab)) == 1){
    mitab_t = unique(copy(mitab[,.(IDs_interactor_A, IDs_interactor_B)]))
    degrees = mitab_t[, .(ID = c(IDs_interactor_A, IDs_interactor_B))][, .N, by = ID]
  } else stop("IDs_interactor_A and IDs_interactor_B columns not supplied")
  if(order_by_degree) setorder(degrees, N, ID)
  return(degrees)
}
