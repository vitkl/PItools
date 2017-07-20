##' find the degree for each protein given molecular interaction data containing unique pair_id
##' @name edgelist2degree
##' @author Vitalii Kleshchevnikov
##' @details unique pair_id is the unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs \code{\link{cleanMITAB}}
##' @param mitab data.table containing molecular interaction data, including pair_id
##' @param sep how are individual interactors in the pair_id reparated? default: "\\|"
##' @return sorted data.table containing degree for each of the nodes in the input, details: \code{mitab}
##' @import data.table
##' @export edgelist2degree
##' @examples
##' edgelist = data.table(pair_id = paste0(rep(c("P", "Q", "R"), each = 30),
##'                       "|",
##'                       sample(rep(c("Z", "X", "F"), each = 30)), 1:90))
##' edgelist2degree(edgelist)
edgelist2degree = function(mitab, sep = "\\|"){
  if("pair_id" %in% colnames(mitab)){
    mitab_t = unique(copy(mitab[,.(pair_id)]))
    mitab_t[, c("ida", "idb") := tstrsplit(pair_id, sep)]
    degrees = mitab_t[, .(c(ida, idb))][, .N, by = V1]
    setnames(degrees, colnames(degrees), c("ID", "N"))
  } else stop("pair_id column not supplied")
  setorder(degrees, N, ID)
  return(degrees)
}
