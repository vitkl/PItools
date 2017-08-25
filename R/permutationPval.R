##' calculate permutation pval for association between X and Z if both are connected to Y
##' @name permutationPval
##' @param interactions2permute formula specifying columns that contain X-Y interactions that will be permuted
##' @param associations2test formula specifying columns that contain X-Z correspondence: empirical p-values for the association of X with Z will be tested, NOTE: Y-Z interactions are assumed
##' @param node_attr formula or list of formulas specifying columns that contain attributes of \code{X}, \code{Y} or \code{Z} or their combination (\code{X} ~ \code{degree}, \code{X} + \code{Z} ~ \code{pvalue})
##' @param data data.table containing interaction data and attributes
##' @param statistic formula that specifies how to calculate statistic using attibutes from \code{node_attr} by node \code{X} in \code{interactions2permute} and \code{Z} in \code{associations2test}, details: \code{\link[MItools]{permutationPvalHelper}}. In data.table synthax: \code{DT[, (observed/permuted)statistic := eval(right-hand-side of formula), by = .(eval(column names in the left-hand-side of formula))]}
##' @param select_nodes formula or list of formulas specifying which nodes of specific node type to select before permutation based on condition (\code{X} ~ \code{degree} > \code{10})
##' @param N number of times to run permutation of PPI network
##' @param cores specify how many cores to use for parallel processing, default (NULL) is to detect all cores on the machine and use all minus one. When using LSF cluster you must specify the number of cores to use because \code{\link[parallel]{detectCores}} doen't know how much cores you have requested from LSF (with bsub -n) and detects all cores on the actual physical node.
##' @param seed seed for RNG for reproducible sampling
##' @param also_permuteYZ logical, permute Y-Z interactions in addition to X-Y (specified in interactions2permute) ?
##' @import data.table
##' @import qvalue
##' @import BiocGenerics
##' @author Vitalii Kleshchevnikov
##' @export permutationPval
permutationPval = function(interactions2permute = nodeX ~ nodeY, associations2test = nodeX ~ nodeZ, node_attr = NULL, data, statistic, select_nodes = NULL, N = 1000, cores = NULL, seed = NULL, also_permuteYZ = F){
  ########################################################################################################################
  # if data is not data.table or is not coerce-able to data.table: stop
  if(!is.data.table(data)) if(is.data.frame(data)) data = as.data.table(data) else if(is.matrix(data)) data = as.data.table(data) else stop("data is provided but is not data.table, data.frame or matrix")
  # if character columns contain empty character "" => substutute that for NA
  col_class_character = sapply(data, function(data_col) class(data_col)) == "character"
  col_character = names(col_class_character)[col_class_character]
  for (colname in col_character) {
    data[eval(as.formula(paste0("~ ", colname," == \"\""))[[2]]), c(colname) := NA]
  }
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
  interactionsXY_cols = c(nodes$nodeX, nodes$nodeY)
  interactionsYZ_cols = c(nodes$nodeY, nodes$nodeZ)
  associations_cols = c(nodes$nodeX, nodes$nodeZ)
  cols = list(interactionsXY_cols = interactionsXY_cols, interactionsYZ_cols = interactionsYZ_cols, associations_cols = associations_cols)
  #extract node attributes
  cols = node_attr2colnames(node_attr, cols, nodes)
  # create tables as specified
  interactionsXY = unique(data[, cols$interactionsXY_cols, with = F])
  interactionsYZ = unique(data[, cols$interactionsYZ_cols, with = F])
  associations = unique(data[, cols$associations_cols, with = F])
  data_list = list(interactionsXY = interactionsXY, interactionsYZ = interactionsYZ, associations = associations)
  ########################################################################################################################
  # check if select_nodes asks to select nodes based on their attributes as specified in node_attr or (by the node name)


  # filter tables by node attribute if select_nodes is provided
  data_list = filterByList(data_list, select_nodes, cols, nodes)
  ########################################################################################################################

  # extract "by columns" and how to calculate statistic from formula provided in statistic argument as class "call"
  by_cols = as.formula(paste0("~ .(",paste0(all.vars(statistic[[2]]), collapse = ","),")"))[[2]]
  exprs = statistic[[3]]
  exprs_vars = all.vars(exprs)
  # check if any X-Z parameteters are needed to calculate statistic
  if(sum(exprs_vars %in% cols$associations_cols) >= 1) includeAssociations = T else includeAssociations = F

  # calculate observed statistic
  data_list = calcObservedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations)

  # create a table to be used for permuting interactions
  data_list$permuted_interactionsXY = data_list$interactionsXY
  # if both XY and YZ interactions are to be permuted create a table to be used for permuting YZ interactions
  if(also_permuteYZ) data_list$permuted_interactionsYZ = data_list$interactionsYZ

  # set up parallel processing
  # create cluster
  if(is.null(cores)) cores = detectCores()-1
  cl <- makeCluster(cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(MItools); library(data.table); library(BiocGenerics)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("data_list", "by_cols", "exprs", "nodes", "nodes_call", "includeAssociations", "also_permuteYZ"), envir=environment())
  # set seed
  clusterSetRNGStream(cl, iseed = seed)

  # splitting computation into the outer and inner replicate helps to save memory by decreasing the total size of the result
  if(N%%10 == 0) {
    if(N < 10000) outer_N = N/10
    if(N < 10000) inner_N = 10
    if(N > 10000 & N%%100 == 0) outer_N = N/100 else outer_N = N/10
    if(N > 10000 & N%%100 == 0) inner_N = 100 else inner_N = 10
  } else {outer_N = N; inner_N = 1}
  clusterExport(cl, c("inner_N"), envir=environment())

  #perform permutations - returns counts when observed statistic is lower than permuted (higher_counts), how many X-Z pair have non missing values and were used in calculation (not_missing), + names of X and Z
  # outer replicate
  temp = parReplicate(cl = cl, n = outer_N, expr = {
    # inner replicate
    temp_inner = replicate(n = inner_N, expr = {
      # calculate statistic using permuted network
      data_list = MItools:::calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ)
      # count how many times observed is lower than permuted giving us the empirical probability of observing value as high or higher by chance
      data_list_temp = MItools:::observedVSpermuted(data_list, nodes_call, nodes)
    })
    # end of inner replicate
    # aggregate attributes across permutations
    temp2_inner = MItools:::aggregatePermutations(temp_inner, nodes, nodes_call)
  }, simplify = TRUE, USE.NAMES = TRUE)
  # end of outer replicate
  stopCluster(cl)

  # aggregate attributes across permutations
  temp2 = aggregatePermutations(temp, nodes, nodes_call)

  # calculate P-value
  temp2[, p.value := higher_counts / not_missing]
  # merge p-value result to the original "data" data.table
  data_with_pval = temp2[data, on = c(nodes$nodeX, nodes$nodeZ)]
  # add observed statistic to the original "data" data.table
  data_list$observed = data_list$observed[,.(nodes$nodeX, nodes$nodeZ, "observed_statistic", "YmissingZ_perX")]
  setnames(data_list$observed, colnames(data_list$observed), c(nodes$nodeX, nodes$nodeZ, "observed_statistic", "YmissingZ_perX"))
  data_with_pval = data_list$observed[data, on = c(nodes$nodeX, nodes$nodeZ)]


  return(list(input = match.call(), nodes = nodes, data_with_pval = data_with_pval))

}

#data = fread("../viral_project/processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
#permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human,
#                associations2test = IDs_interactor_viral ~ IDs_domain_human,
#                node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
#                                 IDs_domain_human ~ domain_count,
#                                 IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
#                data = data,
#                statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
#                select_nodes = IDs_domain_human ~ domain_count > 16,
#                N = 10,
#                cores = NULL, seed = 1)
# fisher.test(matrix(c(data$domain_count[1], data$N_prot[1] - data$domain_count[1], data$domain_count_per_IDs_interactor_viral[1], data$IDs_interactor_viral_degree[1] - data$domain_count_per_IDs_interactor_viral[1]),2,2), alternative = "greater", conf.int = F)$estimate
# statistic = IDs_interactor_viral + IDs_domain_human ~ fisher.test(matrix(c(domain_count[1], N_prot[1] - domain_count[1], domain_count_per_IDs_interactor_viral[1], IDs_interactor_viral_degree[1] - domain_count_per_IDs_interactor_viral[1]),2,2), alternative = "greater", conf.int = F)$estimate #p.value
