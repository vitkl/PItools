##' calculate permutation pval for association between X and Z if both are connected to Y
##' @name permutationPval
##' @param interactions2permute
##' @param associations2test
##' @param node_attr formula or list of formulas specifying columns that contain attributes of \code{X}, \code{Y} or \code{Z} or their combination (\code{X} ~ \code{degree}, \code{X} + \code{Z} ~ \code{pvalue})
##' @param data
##' @param statistic function that calculates statistic which takes columns specified in \code{interactions2permute}, \code{associations2test} or \code{node_attr}
##' @param select_nodes formula or list of formulas specifying which nodes of specific node type to select before permutation based on condition (\code{X} ~ \code{degree} > \code{10})
##' @param N number of times to run permutation of PPI network
##' @param cores specify how many cores to use for parallel processing, default (NULL) is to detect all cores on the machine and use all minus one. When using LSF cluster you must specify the number of cores to use because \code{\link[BiocGenerics]{detectCores}} doen't know how much cores you have requested from LSF (with bsub -n) and detects all cores on the actual physical node.
##' @param seed seed for RNG for reproducible sampling
##' @import data.table
##' @import qvalue
##' @author Vitalii Kleshchevnikov
##' @export permutationPval
permutationPval = function(interactions2permute = nodeX ~ nodeY, associations2test = nodeX ~ nodeZ, node_attr = NULL, data, statistic, select_nodes = NULL, N = 1000, cores = NULL, seed = 1, include_missing_Z_as_zero = F){
  ########################################################################################################################
  # if data is not data.table or is not coerce-able to data.table: stop
  if(!is.data.table(data)) if(is.data.frame(data)) data = as.data.table(data) else if(is.matrix(data)) data = as.data.table(data) else stop("data is provided but is not data.table, data.frame or matrix")
  ########################################################################################################################
  # extract nodes from formulas: interactions2permute and associations2test
  nodes = formula2nodes(interactions2permute, associations2test)
  # a way to convert character column name into code executable by data.table calls
  nodes_call = list()
  nodes_call$nodeX = as.formula(paste0("~ ", nodes$nodeX))[[2]]
  nodes_call$nodeY = as.formula(paste0("~ ", nodes$nodeY))[[2]]
  nodes_call$nodeZ = as.formula(paste0("~ ", nodes$nodeZ))[[2]]
  ########################################################################################################################
  # find columns in data that should be separated into individual data.table-s specified in interactions2permute, associations2test and node_attr
  interactionsXY_cols = c(nodeX, nodeY)
  interactionsYZ_cols = c(nodeY, nodeZ)
  associations_cols = c(nodeX, nodeZ)
  cols = list(interactionsXY_cols = interactionsXY_cols, interactionsYZ_cols = interactionsYZ_cols, associations_cols = associations_cols)
  #extract node attributes
  cols = node_attr2colnames(node_attr, cols, nodes)
  # create tables as specified
  interactionsXY = unique(data[, cols$interactionsXY_cols, with = F])
  interactionsYZ = unique(data[, cols$interactionsYZ_cols, with = F])
  associations = unique(data[, cols$associations_cols, with = F])
  data_list = list(interactionsXY = interactionsXY, interactionsYZ = interactionsYZ, associations = associations)
  ########################################################################################################################
  # filter tables by node attribute if select_nodes is provided
  # extract nodes to filter and apply conditions
  if(is.list(select_nodes)){
    # how many formulas in a list?
    N_attr = length(select_nodes)
    list_names = character(N_attr)
    for(i in 1:N_attr){ # for each formula extract elements
      form_temp = select_nodes[[i]]
      list_names[i] = as.character(as.expression(form_temp[[2]]))

      data_list = filterByFormula(data_list, form_temp, cols, nodes)
    }
    # give list elements names
    names(select_nodes) = list_names
  } else if(is.formula(select_nodes)){ # extract element from single formula
    data_list = filterByFormula(data_list, select_nodes, cols, nodes)
  } else if(is.null(node_attr)) NULL else stop("node_attr is provided but is neither a list nor a formula")
  ########################################################################################################################

  # extract "by columns" and how to calculate statistic from formula provided in statistic argument as class "call"
  by_cols = as.formula(paste0("~ .(",paste0(all.vars(statistic[[2]]), collapse = ","),")"))[[2]]
  exprs = statistic[[3]]

  # calculate observed statistic
  data_list = calcObservedStatistic(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero)

  data_list = calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero)
}

data[eval(as.expression((~ domain_count == 4)[[2]])),]

data = fread("../viral_project/processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
                associations2test = IDs_interactor_viral ~ IDs_domain_human,
                node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                 IDs_domain_human ~ domain_count,
                                 IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                data = data,
                statistic = function(x) x^2,
                select_nodes = IDs_domain_human ~ domain_count > 16,
                N = 1000,
                cores = NULL, seed = 1)
