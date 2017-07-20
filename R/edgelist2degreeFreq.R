##' Find the degree distribution given edgelist (pair_id) and a list of proteins which degrees matter
##' @name edgelist2degreeFreq
##' @author Vitalii Kleshchevnikov
##' @details unique pair_id is the unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs \code{\link{cleanMITAB}}
##' @param mitab data.table containing molecular interaction data, including pair_id
##' @param prots list of relevant proteins, degree distribution is based on the prevalence of specific degrees in this list, all proteins in \code{mitab} by default
##' @param sep ID separator in pair_id, passed to \code{\link{edgelist2degree}}
##' @import data.table
##' @export edgelist2degreeFreq
##' @examples
##' set.seed(1)
##' random = randomInteractome(n_prot = 200, degree_dist = NULL, taxid = "9606", database = "imex", protein_only = TRUE)
##' edgelist2degreeFreq
edgelist2degreeFreq = function(mitab, prots = NULL, sep = "\\|"){
  if(!(is.data.table(mitab) & "pair_id" %in% colnames(mitab))) stop("mitab table (edge list) is not data.table or doesn't contain pair_id column")
  # use all proteins in mitab if no protein list supplied
  if(is.null(prots)) prots = mitab[, unique(unlist(strsplit(pair_id, "\\|")))]

  mitab = copy(mitab)
  # compute degree
  degree = edgelist2degree(mitab, sep)[ID %in% prots,]
  # compute degree count
  degree_dist = unique(degree[, count := .N, by = N][, .(N,count)])
  # compute degree frequency
  degree_dist[, degree_freq := count / sum(count)]
  return(degree_dist[,.(N, degree_freq)])
}
