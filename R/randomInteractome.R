##' Retrieve molecular interactions for the random set of proteins (of a particular taxon)
##' @name randomInteractome
##' @author Vitalii Kleshchevnikov
##' @param MITABdata data.table containing pre-loaded molecular interaction data as returned by \code{\link{queryPSICQUICrlib}}, usefull for taking multiple samples, the default in NULL
##' @param degree_data data.table containing pre-calculated (using \code{\link{edgelist2degree}}) degree for each node in MITABdata, usefull for taking multiple samples, the default in NULL
##' @param n_prot integer (1L), the number of proteins for which to retrieve the random set of interactions
##' @param degree_dist data.table, specifies the degree frequency () for each degree (N) to produce the network with the specific degree distribution, if set to NULL (default) the degree distribution will correspond to that of \code{taxid} interactome
##' @return list of two elements: 1. interactome - data.table containing molecular interaction data in either of these two formats:
##' @return if \code{clean} is TRUE: contains columns as described in \code{\link{cleanMITAB}};
##' @return if \code{clean} is FALSE: contains a standard set of columns for MITAB2.5 or MITAB2.7 depending on \code{format};
##' @return 2. seed - character vector containing IDs of proteins used as a seed to retrieve molecular interactions
##' @details Random network can be specified to have specific degree distribution. If the (\code{degree} parameter is set \code{taxid} proteins will be split by degree and from each degree group a sample of the size specified by how many times specific degree number is repeated in \code{degree} will be taken.
##' @details If the degree distribution is not specified a sample of \code{n_prot} is taken from all proteins which have interaction data available in the \code{database} for \code{taxid}. In this case, the degree distribution of the resulting set of proteins will be similar to the degree distribution in the interactome of \code{taxid} in \code{database}.
##' @details \code{randomInteractome} retrieves molecular interactions using \code{\link{fullInteractome}}
##' @import data.table
##' @export randomInteractome
##' @examples
##' # retrive the interactome using PSICQIUC servise (or by reading local copy) from IMEx databases for a list of 200 random human (9606) proteins, not specifying their degree distribution
##' set.seed(1)
##' random = randomInteractome(n_prot = 200, degree_dist = NULL, taxid = "9606", database = "imex", protein_only = TRUE)
##' # retrive the interactome from MITABdata for a list of 200 random human (9606) proteins, not specifying their degree_dist distribution
##' full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE)
##' set.seed(1)
##' random = randomInteractome(MITABdata = full, n_prot = 200, degree_dist = NULL)
##' # retrive the interactome using PSICQIUC servise (or by reading local copy from a specified directory) from IMEx databases for a list of 200 random human (9606) proteins, not specifying their degree_dist distribution
##' set.seed(1)
##' random = randomInteractome(n_prot = 200, degree_dist = NULL, taxid = "9606", database = "imex", protein_only = TRUE, directory = "./data/")
randomInteractome = function(MITABdata = NULL, degree_data = NULL, n_prot, degree_dist = NULL, ...){
  # if MITABdata is NULL retrive the cleaned full Interactome
  if(is.null(MITABdata)) full_interactome_clean = fullInteractome(..., format = "tab25", clean = TRUE)
  # if MITABdata is supplied use that data
  if(!is.null(MITABdata)) full_interactome_clean = MITABdata

  # check if the data has necessary columns:
  if(mean(c("pair_id") %in% colnames(full_interactome_clean)) != 1) stop("MITABdata is in the wrong format: no pair_id column")

  # get interactors
  interactors = full_interactome_clean[, unique(unlist(strsplit(pair_id, "\\|")))]

  # if the degree_dist distribution is not specified - sample n_prot of interactors
  if(is.null(degree_dist)){
    random_interactors = sample(interactors, n_prot)
  }

  # if the degree_dist distribution is specified - calculate degree_dist of each interactor (the number of interacting partners per interactor)
  if(!is.null(degree_dist)){
    if(mean(c("N", "degree_freq") %in% colnames(degree_dist)) != 1) stop("degree_dist does not contain \"N\" and \"degree_freq\" columns")
    # calculate degree_data data if not provided
    if(is.null(degree_data)) degree_data = edgelist2degree(full_interactome_clean[,.(pair_id)])
    # stop if degree_data doesn't match MITABdata
    if(!is.null(degree_data) & mean(degree_data$ID %in% interactors) != 1) stop("degree_data doesn't contain degree_dist information for all nodes in MITABdata")
    # attach each interactor degree_dist for the sample
    degree_data = degree_data[degree_dist, on = "N"]
    random_interactors = sample(x = degree_data[degree_freq != 0, ID],
                                size = n_prot,
                                prob = degree_data[degree_freq != 0, degree_freq])
  }

  # retrive interactions for a set of random proteins
  return(list(interactome = interactors2interactions(full_interactome_clean, random_interactors), seed_proteins = random_interactors))

}
