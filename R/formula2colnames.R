##' functions that are called inside permutationPval
##' @name permutationPvalHelper
##' @rdname permutationPvalHelper
##' @export permutationPvalHelper
permutationPvalHelper = NULL
##' extract and match attribute column names to data.tables that contain nodes which those attributes describe
##' @name formula2colnames
##' @rdname permutationPvalHelper
##' @param inheritParams permutationPval node_attr interactions2permute associations2test select_nodes also_permuteYZ
##' @param cols list of column names
##' @param nodes list of node names (column names for columns that contain node names)
##' @param nodes_call list of node names, each class call (to be evaluated within data.table expressions DT[,eval(nodes_call$nodeX)])
##' @param data_list list that contains data.tables necessary for these functions to work
##' @param by_cols modified left-hand side of \code{statistic} (\code{\link[MItools]{permutationPval}}), class call
##' @param exprs right-hand side of \code{statistic} (\code{\link[MItools]{permutationPval}}), class call, \code{exprs} is evaluated by \code{by_cols}; in data.table synthax: DT[, statistic := eval(exprs), by = .(eval(by_cols))]
##' @param includeAssociations logical, if calculating statistic requires columns that contain attribute of X-Z pair associations table will be merged
##' @details extract nodes
##' @details extract attributes
##' @details attach attributes to a table if all nodes are present in that table
##' @details and in case of associations table both varibles in a formula should be nodeX and nodeZ
##' @return list of column names for XY, YZ and association tables
##' @import data.table
##' @import BiocGenerics
##' @author Vitalii Kleshchevnikov
##' @usage cols = formula2colnames(node_attr, cols, nodes)
##' nodes_vs_attributes = formula2colnames(node_attr, nodes_vs_attributes = T)
##' data_list = calcObservedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations)
##' data_list = calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ)
##' res = observedVSpermuted(data_list, nodes_call, nodes)
##' temp2 = aggregatePermutations(temp, nodes, nodes_call)
formula2colnames = function(node_attr, cols, nodes, nodes_vs_attributes = F){
  # if not formula is provided
  if(!is.formula(node_attr)) stop(paste("node_attr is provided but is not a formula: class - ", class(node_attr), "; content - ",paste(paste0("[[",1:length(node_attr),"]]"), node_attr, collapse = " ")))

  # extract nodes
  node_temp = all.vars(node_attr[[2]])
  # extract attributes
  attribute_temp = all.vars(node_attr[[3]])
  if(!nodes_vs_attributes){
    # attach attributes to a table if all nodes are present in that table
    if(mean(node_temp %in% c(nodes$nodeX, nodes$nodeY)) == 1) cols$interactionsXY_cols = c(cols$interactionsXY_cols, attribute_temp)
    if(mean(node_temp %in% c(nodes$nodeY, nodes$nodeZ)) == 1) cols$interactionsYZ_cols = c(cols$interactionsYZ_cols, attribute_temp)
    # and in case of associations table both varibles in a formula should be nodeX and nodeZ
    if(mean(c(nodes$nodeX, nodes$nodeZ) %in% node_temp) == 1) cols$associations_cols = c(cols$associations_cols, attribute_temp)
    return(cols)
  }
  if(nodes_vs_attributes) return(nodes = node_temp, attributes = attribute_temp)
}

##' @name node_attr2colnames
##' @rdname permutationPvalHelper
##' @return the same as \code{formula2colnames} but when node_attr2colnames is a list of formulas
##' @usage cols = node_attr2colnames(node_attr, cols, nodes)
node_attr2colnames = function(node_attr, cols, nodes) {

  if(is.list(node_attr)){
    # how many formulas in a list?
    N_attr = length(node_attr)
    list_names = character(N_attr)
    for(i in 1:N_attr){ # for each formula extract elements
      form_temp = node_attr[[i]]
      list_names[i] = as.character(as.expression(form_temp[[2]]))
      cols = formula2colnames(node_attr = form_temp, cols, nodes)
    }
    # give names to the elements of the list
    names(node_attr) = list_names
  } else if(is.formula(node_attr)){ # extract element from single formula
    cols = formula2colnames(node_attr = node_attr, cols, nodes)
  } else if(is.null(node_attr)) NULL else
    stop("node_attr is provided but is neither a list nor a formula")

  return(cols)
}

##' @name formula2nodes
##' @rdname permutationPvalHelper
##' @details \code{formula2nodes} reorders node column names so that regardless of their order in formula we later permute XY interactions and test for XZ associations
##' @return list of 3 node column names extracted from formulas supplied through \code{interactions2permute}, \code{associations2test}
##' @usage nodes = formula2nodes(interactions2permute, associations2test)
formula2nodes = function(interactions2permute, associations2test){
  # if not formulas are provided
  if(!is.formula(interactions2permute)) stop("interactions2permute is provided but is not a formula")
  if(!is.formula(associations2test)) stop("associations2test is provided but is not a formula")

  # extract nodes to permute from formula ([[2]] extracts the first element of the formula, [[3]] extracts the second)
  nodeX = all.vars(interactions2permute[[2]])
  nodeY = all.vars(interactions2permute[[3]])
  # extract nodes for which to test association from formula
  nodeXtest = all.vars(associations2test[[2]])
  nodeYtest = all.vars(associations2test[[3]])
  # generate standard set of variable names: find the name of the thing that is asked to be associated with nodeX
  if(nodeX == nodeXtest) nodeZ = nodeYtest
  if(nodeX == nodeYtest) nodeZ = nodeXtest
  # generate standard set of variable names: find the name of the thing that is asked to be associated with nodeY
  if(nodeY == nodeXtest) nodeZ = nodeYtest
  if(nodeY == nodeYtest) nodeZ = nodeXtest
  # if formula suggest to find association between nodeY and some other column - reorder
  if(nodeY == nodeXtest | nodeY == nodeYtest) {
    nodeX_temp = nodeY
    nodeY = nodeX
    nodeX = nodeX_temp
  }
  # check if all 3 are different
  if(nodeX == nodeY | nodeX == nodeZ | nodeY == nodeZ) stop(paste0("formulas interactions2permute and associations2test should provide different input, such as: \n interactions2permute = nodeX ~ nodeY, associations2test = nodeX ~ nodeZ"))

  nodes = list(nodeX = nodeX, nodeY = nodeY, nodeZ = nodeZ)
  return(nodes)
}

##' @name filterByFormula
##' @rdname permutationPvalHelper
##' @details filters one or 2 of the XY, YZ or association tables based on a condition defining which nodes to keep (only one node type: X, Y or Z)
##' @return data_list in which interactionsXY, interactionsYZ, associations data.table-s were filtered by formula
##' @usage data_list = filterByFormula(data_list, select_nodes, cols, nodes)
filterByFormula = function(data_list, select_nodes, cols, nodes){
  # if not formula is provided
  if(!is.formula(select_nodes)) stop("select_nodes is provided but is not a formula")

  # extract node column
  node_temp = all.vars(select_nodes[[2]])
  # extract node column (class name)
  node_ftemp = select_nodes[[2]]
  # extract condition (class call)
  condition_ftemp = select_nodes[[3]]
  # extract condition attribute names
  attribute_temp = all.vars(select_nodes[[3]])
  if(length(node_temp) > 1) stop("selecting multiple nodes by single condition is not implemented yet, however, you can cheat by supplying the same condition twice for two different node names")

  # find which nodes to keep
  # finding which data.table contains attribute necessary for filtering AND has node type that we want to filter ->
  # -> then creates a vector of nodeIDs of specific node type by evaluating the expression on the right hand side of the formula
  if(mean(attribute_temp %in% cols$interactionsXY_cols) == 1 & mean(node_temp %in% c(nodes$nodeX, nodes$nodeY)) == 1) node2keep = data_list$interactionsXY[eval(condition_ftemp), eval(node_ftemp)]
  if(mean(attribute_temp %in% cols$interactionsYZ_cols) == 1 & mean(node_temp %in% c(nodes$nodeY, nodes$nodeZ)) == 1) node2keep = data_list$interactionsYZ[eval(condition_ftemp), eval(node_ftemp)]
  if(mean(attribute_temp %in% cols$associations_cols) == 1 & mean(node_temp %in% c(nodes$nodeX, nodes$nodeZ)) == 1) node2keep = data_list$associations[eval(condition_ftemp), eval(node_ftemp)]

  # filter data.table using condition if the node (nodes) are present in that table
  if(mean(node_temp %in% c(nodes$nodeX, nodes$nodeY)) == 1) data_list$interactionsXY = data_list$interactionsXY[eval(node_ftemp) %in% node2keep,]
  if(mean(node_temp %in% c(nodes$nodeY, nodes$nodeZ)) == 1) data_list$interactionsYZ = data_list$interactionsYZ[eval(node_ftemp) %in% node2keep,]
  if(mean(node_temp %in% c(nodes$nodeX, nodes$nodeZ)) == 1) data_list$associations = data_list$associations[eval(node_ftemp) %in% node2keep,]
  return(data_list)
}

calcObservedStatistic = function(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations){

  # merge without removing X that don't have a match in Z, interactionsXY on the inside of "[" (in i position) means keep all interactionsXY, discard non-matching interactionsYZ
  data_list$observed = unique(data_list$interactionsXY[data_list$interactionsYZ, on = nodes$nodeY, allow.cartesian = T, nomatch = 0])
  # if calculating statistic requires some parameters of both X and Z -> merge associations table containing necessary data
  if(includeAssociations) data_list$observed = unique(data_list$associations[data_list$observed, on = c(nodes$nodeX, nodes$nodeZ), allow.cartesian = T, nomatch = 0])

  # record in how many cases statistic is missing per X
  data_list$observed[, YmissingZ_perX := sum(is.na(eval(nodes_call$nodeZ))),
                     by = eval(nodes_call$nodeX)]
  #NAs are removed otherwise statistic cannot be calculated
  data_list$observed = data_list$observed[!is.na(eval(nodes_call$nodeZ)),]

  # calculate observed statistic
  data_list$observed[, observed_statistic := eval(exprs), by = eval(by_cols)]

  return(data_list)
}

calcPermutedStatistic = function(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ){

  # permute XY network (note: [, class character (1L) := eval(class call)])
  data_list$permuted_interactionsXY[, nodes$nodeY := sample(eval(nodes_call$nodeY))]
  # if also_permuteYZ is TRUE permute YZ network
  if(also_permuteYZ) data_list$permuted_interactionsYZ[, nodes$nodeY := sample(eval(nodes_call$nodeY))]

  # merge without removing X that don't have a match in Z, interactionsXY on the inside of "[" (in i position) means keep all interactionsXY, discard non-matching interactionsYZ
  # if also_permuteYZ is FALSE merge permuted XY to observed YZ
  if(!also_permuteYZ) data_list$permuted = unique(data_list$permuted_interactionsXY[data_list$interactionsYZ, on = nodes$nodeY, allow.cartesian = T, nomatch = 0])
  # if also_permuteYZ is FALSE merge permuted XY to permuted YZ
  if(also_permuteYZ) data_list$permuted = unique(data_list$permuted_interactionsXY[data_list$permuted_interactionsYZ, on = nodes$nodeY, allow.cartesian = T, nomatch = 0])
  # if calculating statistic requires some parameters of both X and Z -> merge associations table containing necessary data
  if(includeAssociations) data_list$permuted = unique(data_list$associations[data_list$permuted, on = c(nodes$nodeX, nodes$nodeZ), allow.cartesian = T, nomatch = 0])

  # record in how many cases statistic is missing per X: we don't care in permuted cases, so I don't track this
  # data_list$permuted[, YmissingZ_perX := sum(is.na(eval(nodes_call$nodeZ))),
  #                    by = eval(nodes_call$nodeX)]
  #NAs are removed otherwise statistic cannot be calculated
  data_list$permuted = data_list$permuted[!is.na(eval(nodes_call$nodeZ)),]

  # calculate permuted statistic
  data_list$permuted[, permuted_statistic := eval(exprs), by = eval(by_cols)]

  # since this statistic is being calculated to generate background distribution for nodeX - we keep only the node X and the statistic

  # collapse (unique(), to remove X-Z pairs duplicated because of many Y) counts by nodeX AND nodeZ, then remove nodeZ: otherwise identical values of counts will get collapsed into one (for example, when two nodeZ have the same frequency among interactors of one nodeX)
  net_notNA = unique(data_list$permuted[!is.na(eval(nodes_call$nodeZ)), .(eval(nodes_call$nodeX), eval(nodes_call$nodeZ), permuted_statistic)])
  setnames(net_notNA, colnames(net_notNA), c(nodes$nodeX, nodes$nodeZ,"permuted_statistic"))
  data_list$permuted = net_notNA[, nodes$nodeZ := NULL]

  return(data_list)
}

observedVSpermuted = function(data_list, nodes_call, nodes){
  # merge keeping all observed in i and discarding all permuted in x
  result = data_list$permuted[data_list$observed, on = nodes$nodeX, allow.cartesian=TRUE]
  # calculate in how many cases permuted statistic is at least as large as observed statistic
  result_pval = result[, sum(observed_statistic <= permuted_statistic, na.rm = T), by = .(eval(nodes_call$nodeX), eval(nodes_call$nodeZ))]
  setnames(result_pval, colnames(result_pval), c(nodes$nodeX, nodes$nodeZ, "higher_counts"))
  setorderv(result_pval, cols = c(nodes$nodeX, nodes$nodeZ, "higher_counts"))

  # calculate in how many cases there is a match between X and Z => count of (permutations * number of comparison in each permutation) per each X-Z pair
  result_not_missing = result[, sum(!(is.na(observed_statistic) | is.na(permuted_statistic))), by = .(eval(nodes_call$nodeX), eval(nodes_call$nodeZ))]
  setnames(result_not_missing, colnames(result_not_missing), c(nodes$nodeX, nodes$nodeZ, "not_missing"))
  setorderv(result_not_missing, cols = c(nodes$nodeX, nodes$nodeZ, "not_missing"))

  #merge these higher_counts and not_missing
  temp = result_pval[result_not_missing, on = c(nodes$nodeX, nodes$nodeZ)][!is.na(eval(nodes_call$nodeZ)),]

  return(temp[, c("higher_counts", "not_missing", nodes$nodeX, nodes$nodeZ), with = F])
}

aggregatePermutations = function(temp, nodes, nodes_call){
  # aggregate attributes across permutations
  temp2 = data.table(higher_counts = 0, not_missing = 0, temp[nodes$nodeX,1][[1]], temp[nodes$nodeZ,1][[1]])
  setnames(temp2, colnames(temp2), c("higher_counts", "not_missing", nodes$nodeX, nodes$nodeZ))

  for(i in 1:ncol(temp)){
    if(!(all.equal(temp2[, eval(nodes_call$nodeX)], temp[nodes$nodeX, i][[1]]) &
         all.equal(temp2[, eval(nodes_call$nodeZ)], temp[nodes$nodeZ, i][[1]]))) stop("X and/or Z ids mismatch, you have stumbled upon a bug in permutationPval code") else
           temp2[, higher_counts := higher_counts + temp["higher_counts", i][[1]]]
  }
  for(i in 1:ncol(temp)){
    temp2[, not_missing := not_missing + temp["not_missing",i][[1]]]
  }
  return(temp2)
}
