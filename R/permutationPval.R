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
##' @param cluster_type Type of the cluster to create when using R parrallel (\code{\link[parallel]{clusterApply}}). Type "FORK" means cluster nodes share objects in memory. Details: (\code{\link[parallel]{makeCluster}})
##' @param seed seed for RNG for reproducible sampling
##' @param also_permuteYZ logical, permute Y-Z interactions in addition to X-Y (specified in interactions2permute) ?
##' @param clustermq if TRUE uses clustermq job scheduling (\code{\link[clustermq]{Q}}) instead of local parallelisation (\code{\link[MItools]{parReplicate}})  = F,  = 4000,  = 100, clustermq_template = list(), split_comp_inner_N
##' @param clustermq_mem memory in MB to allocate for each job (ignored unless clustermq == TRUE)
##' @param clustermq_jobs maximal number of computing cluster jobs to use (ignored unless clustermq == TRUE)
##' @param clustermq_template Add specific arguments to computing cluster job submission call. Not needed in most cases. Details: \code{\link[clustermq]{Q}} (ignored unless clustermq == TRUE)
##' @param split_comp_inner_N parallel evaluation of permutations is split into the outer and inner replicate calls helps to save memory by decreasing the total size of the result. This argument let's you manually specify the number of inner replicate calls. This has to be optimised for data size when using clustermq (ignored unless clustermq == TRUE)
##' @param clustermq_fail_on_error If TRUE clustermq will fail if one of the jobs returns an error. Details: \code{\link[clustermq]{Q}} (ignored unless clustermq == TRUE)
##' @param clustermq_log_worker If TRUE clustermq will save log of worker jobs. Where in is save is determined by clustermq.template
##' @param formula argument for \code{permutationPvalPlot}, formula specifying attribute of which nodes to plot like this: nodeX + nodeZ ~ p.value. The default is to plot p.value histogram for nodeX and nodeZ as specified in the \code{x} object
##' @param x argument for \code{permutationPvalPlot}, output of \code{permutationPval}, class "XYZinteration_XZEmpiricalPval"
##' @param ... argument for \code{permutationPvalPlot}, base R plotting parameters
##' @return object of S3 class "XYZinteration_XZEmpiricalPval" (list), containing \code{permutationPval} function call, standardised node names, and data.table containing the original data but appended with empirical p-value (p.value), observed_statistic, YmissingZ_perX, and higher_counts, not_missing used to calculate p-value
##' @import data.table
##' @import qvalue
##' @import BiocGenerics
##' @author Vitalii Kleshchevnikov
##' @export permutationPval
##' @export plot.XYZinteration_XZEmpiricalPval
##' @export print.XYZinteration_XZEmpiricalPval
##' @usage
##' res = permutationPval(interactions2permute = nodeX ~ nodeY,
##'  associations2test = nodeX ~ nodeZ,
##'  node_attr = NULL, data, statistic,
##'  select_nodes = NULL, N = 1000,
##'  cores = NULL, seed = NULL,
##'  also_permuteYZ = F)
##'
##' # print
##' res
##'
##' # plot p-value distribution (hist)
##' plot(res)
##'
##' # plot the number of Y without Z per X (hist),
##' # formula is used to subset the table before plotting
##' # to avoid plotting single number multiple times
##' plot(res, nodeX ~ YmissingZ_perX)
permutationPval = function(interactions2permute = nodeX ~ nodeY, associations2test = nodeX ~ nodeZ, node_attr = NULL, data, statistic, select_nodes = NULL, N = 1000, cores = NULL, cluster_type = "PSOCK", seed = NULL, also_permuteYZ = F, clustermq = F, clustermq_mem = 4000, clustermq_jobs = 100, clustermq_template = list(), split_comp_inner_N = NULL, clustermq_fail_on_error = TRUE, clustermq_log_worker = FALSE){
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
  if(!is.null(select_nodes)) select_nodes = checkSelectNodesList(select_nodes, node_attr, nodes, cols)

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

  #### parallel processing: splitting computation into the outer and inner replicate helps to save memory by decreasing the total size of the result
  if(N%%10 == 0) {
    if(N < 10000) outer_N = N/10
    if(N < 10000) inner_N = 10
    if(N > 10000 & N%%100 == 0) outer_N = N/100 else outer_N = N/10
    if(N > 10000 & N%%100 == 0) inner_N = 100 else inner_N = 10
  } else {outer_N = N; inner_N = 1}
  if(!is.null(split_comp_inner_N)) {inner_N = split_comp_inner_N; outer_N = round(N/split_comp_inner_N)}

  ########### set up parallel processing
  if(!clustermq) {
    # create cluster
    if(is.null(cores)) cores = detectCores()-1
    cl <- makeCluster(cores, type = cluster_type)
    if(cluster_type != "FORK"){
      # get library support needed to run the code
      clusterEvalQ(cl, {library(MItools); library(data.table); library(BiocGenerics)})
      # put objects in place that might be needed for the code
      clusterExport(cl, c("data_list", "by_cols", "exprs", "nodes", "nodes_call", "includeAssociations", "also_permuteYZ", "inner_N"), envir=environment())
    }
    # set seed
    clusterSetRNGStream(cl, iseed = seed)

    #perform permutations - returns counts when observed statistic is lower than permuted (higher_counts), how many X-Z pair have non missing values and were used in calculation (not_missing), + names of X and Z
    # outer replicate
    temp = parReplicate(cl = cl, n = outer_N, expr = {
      # inner replicate
      temp_inner = replicate(n = inner_N, expr = {
        # calculate statistic using permuted network
        data_list = MItools:::calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ)
        # count how many times observed is lower than permuted giving us the empirical probability of observing value as high or higher by chance
        data_list_temp = MItools:::observedVSpermuted(data_list, nodes_call, nodes)
      }, simplify = FALSE)
      # end of inner replicate
      # aggregate attributes across permutations
      MItools:::aggregatePermutations(temp_inner, nodes, nodes_call)
    }, simplify = FALSE, USE.NAMES = TRUE)
    # end of outer replicate
    stopCluster(cl)
    ########### stop parallel processing
  } else {
    ########### use clustermq
    temp = clustermq::Q(fun = function(i){
      # inner replicate
      temp_inner = replicate(n = inner_N, expr = {
        # calculate statistic using permuted network
        data_list = MItools:::calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ)
        # count how many times observed is lower than permuted giving us the empirical probability of observing value as high or higher by chance
        data_list_temp = MItools:::observedVSpermuted(data_list, nodes_call, nodes)
      }, simplify = FALSE)
      # end of inner replicate
      # aggregate attributes across permutations
      MItools:::aggregatePermutations(temp_inner, nodes, nodes_call)
    }, 1:outer_N,
    export = list(data_list = data_list, by_cols = by_cols,
                  exprs = exprs, nodes = nodes, nodes_call = nodes_call,
                  includeAssociations = includeAssociations,
                  also_permuteYZ = also_permuteYZ, inner_N = inner_N),
    seed = seed,
    memory = clustermq_mem, template = clustermq_template, n_jobs = clustermq_jobs,
    fail_on_error = clustermq_fail_on_error, log_worker = clustermq_log_worker)
  }
  ########### stop clustermq

  # aggregate attributes across permutations
  temp = aggregatePermutations(temp, nodes, nodes_call)

  # calculate P-value
  temp[, p.value := higher_counts / not_missing]
  # merge p-value result to the original "data" data.table
  data_with_pval = temp[data, on = c(nodes$nodeX, nodes$nodeZ)]
  # add observed statistic to the original "data" data.table
  data_list$observed = data_list$observed[, c(nodes$nodeX, nodes$nodeZ, "observed_statistic", "YmissingZ_perX"), with = F]
  data_with_pval = data_list$observed[data_with_pval, on = c(nodes$nodeX, nodes$nodeZ), allow.cartesian = TRUE]

  out = list(input = match.call(), nodes = nodes, data_with_pval = unique(data_with_pval))
  class(out) = "XYZinteration_XZEmpiricalPval"
  return(out)
}

plot.XYZinteration_XZEmpiricalPval = function(x, formula = NULL, main = "", ...){
  if(is.null(formula)) {
    hist(unique(x$data_with_pval[, c(x$nodes$nodeX, x$nodes$nodeZ, "p.value"), with = FALSE])[, p.value],
         breaks = seq(-0.01,1.01,0.01), xlab = "empirical P value", main = main, ...)
  } else if(is.formula(formula)) {
    vars = all.vars(formula[[2]])
    vals = all.vars(formula[[3]])
    hist(unique(x$data_with_pval[, c(vars, vals[1]), with = FALSE])[, eval(formula[[3]])], main = main, ...)
  } else stop("formula argument supplied but is not a formula")
}

print.XYZinteration_XZEmpiricalPval = function(x){
  cat(paste0("\n This object contains the empirical p-value for the association between <", x$nodes$nodeX, "> and <", x$nodes$nodeZ, "> both of which are linked (connected through edges) by <", x$nodes$nodeY, "> \n"))
  cat(paste0("\n Produced by permutationPval function call: \n\n"))
  print(x$input)
  cat("\n $data_with_pval is the data.table containing the original data and appended with empirical p-value (p.value) as well as observed_statistic, YmissingZ_perX, and higher_counts, not_missing used to calculate p-value,\n Details: `?permutationPval` \n\n x$data_with_pval \n")
  print(x$data_with_pval)
}
