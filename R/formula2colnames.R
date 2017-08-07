##' extract and match attribute column names to data.tables that contain nodes which those attributes describe
##' @name formula2colnames
##' @aliases node_attr2colnames
##' @aliases formula2nodes
##' @aliases data_list
##' @param node_attr formula (formula2colnames and node_attr2colnames) or list of formulas (node_attr2colnames)
##' @param cols list of column names
##' @param nodes list of node names (column names for columns that contain node names)
##' @param interactions2permute
##' @param associations2test
##' @param data_list
##' @param select_nodes
##' @param by_cols
##' @param exprs
##' @param include_missing_Z_as_zero
##' @return updated list of column names
##' @details extract nodes
##' @details extract attributes
##' @details attach attributes to a table if all nodes are present in that table
##' @details and in case of associations table both varibles in a formula should be nodeX and nodeZ
##' @import data.table
##' @author Vitalii Kleshchevnikov
##' @usage cols = formula2colnames(node_attr, cols, nodes)
##' cols = node_attr2colnames(node_attr, cols, nodes)
##' nodes = formula2nodes(interactions2permute, associations2test)
##' data_list = filterByFormula(data_list, select_nodes, cols, nodes)
##' calcObservedStatistic(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero)
##' calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero)
formula2colnames = function(node_attr, cols, nodes){
  # if not formula is provided
  if(!is.formula(node_attr)) stop(paste("node_attr is provided but is not a formula: class - ", class(node_attr), "; content - ",paste(paste0("[[",1:length(node_attr),"]]"), node_attr, collapse = " ")))

  # extract nodes
  node_temp = all.vars(node_attr[[2]])
  # extract attributes
  attribute_temp = all.vars(node_attr[[3]])
  # attach attributes to a table if all nodes are present in that table
  if(mean(node_temp %in% c(nodes$nodeX, nodes$nodeY)) == 1) cols$interactionsXY_cols = c(cols$interactionsXY_cols, attribute_temp)
  if(mean(node_temp %in% c(nodes$nodeY, nodes$nodeZ)) == 1) cols$interactionsYZ_cols = c(cols$interactionsYZ_cols, attribute_temp)
  # and in case of associations table both varibles in a formula should be nodeX and nodeZ
  if(mean(c(nodes$nodeX, nodes$nodeZ) %in% node_temp) == 1) cols$associations_cols = c(cols$associations_cols, attribute_temp)
  return(cols)
}

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

calcObservedStatistic = function(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero){
  # merge interactionsXY and interactionsYZ, if not include missing Z for each X as NA (and 0 statistic), then remove any X that don't have a match in Z
  if(!include_missing_Z_as_zero) data_list$observed = unique(data_list$interactionsXY[data_list$interactionsYZ, on = nodes$nodeY, allow.cartesian = T, nomatch = 0])
  # if include missing Z for each X as NA (and 0 statistic) then merge without removing X that don't have a match in Z, interactionsXY on the inside (in i position) means keep all interactionsXY, discard non-matching interactionsYZ
  if(include_missing_Z_as_zero) data_list$observed = unique(data_list$interactionsYZ[data_list$interactionsXY, on = nodes$nodeY, allow.cartesian = T])

  # calculate observed statistic
  data_list$observed[, observed_statistic := eval(exprs), by = eval(by_cols)]

  if(include_missing_Z_as_zero){
    # equate statistic for missing Z cases to 0, record this statistic to a separate column
    data_list$observed[is.na(eval(nodes_call$nodeZ)), missingZ_statistic := observed_statistic]
    data_list$observed[is.na(eval(nodes_call$nodeZ)), observed_statistic := 0]
  }
  return(data_list)
}

calcPermutedStatistic = function(data_list, by_cols, exprs, nodes, nodes_call, include_missing_Z_as_zero){

  # permute XY network (note: [, class character (1L) := eval(class call)])
  data_list$permuted_interactionsXY[, nodes$nodeY := sample(eval(nodes_call$nodeY))]

  # merge permuted_interactionsXY and interactionsYZ, if not include missing Z for each X as NA (and 0 statistic), then remove any X that don't have a match in Z
  if(!include_missing_Z_as_zero) data_list$permuted = unique(data_list$permuted_interactionsXY[data_list$interactionsYZ, on = nodes$nodeY, allow.cartesian = T, nomatch = 0])
  # if include missing Z for each X as NA (and 0 statistic) then merge without removing X that don't have a match in Z, interactionsXY on the inside (in i position) means keep all interactionsXY, discard non-matching interactionsYZ
  if(include_missing_Z_as_zero) data_list$permuted = unique(data_list$interactionsYZ[data_list$interactionsXY, on = nodes$nodeY, allow.cartesian = T])

  # calculate permuted statistic
  data_list$permuted[, permuted_statistic := eval(exprs), by = eval(by_cols)]

  if(include_missing_Z_as_zero){
    # equate statistic for missing Z cases to 0, record this statistic to a separate column
    data_list$permuted[is.na(eval(nodes_call$nodeZ)), missingZ_statistic := permuted_statistic]
    data_list$permuted[is.na(eval(nodes_call$nodeZ)), permuted_statistic := 0]
  }

  # since this statistic is being calculated to generate background distribution for nodeX - we keep only the node X and the statistic

  # separate interactionsXY involving no domains: otherwise these will get collapsed into one per X
  net_NA = data_list$permuted[is.na(eval(nodes_call$nodeZ)), .(eval(nodes_call$nodeX), permuted_statistic)]
  setnames(net_NA, colnames(net_NA), c(nodes$nodeX, "permuted_statistic"))
  # next, collapse (unique()) counts by nodeX AND nodeZ, then remove nodeZ: otherwise identical values of counts will get collapsed into one (for example, when two nodeZ have the same frequency among interactors of one nodeX)
  net_notNA = unique(data_list$permuted[!is.na(eval(nodes_call$nodeZ)), .(eval(nodes_call$nodeX), eval(nodes_call$nodeZ), permuted_statistic)])
  setnames(net_notNA, colnames(net_notNA), c(nodes$nodeX, nodes$nodeZ,"permuted_statistic"))
  net_notNA[, nodes$nodeZ := NULL]

  # rbind no nodeZ and nodeZ present cases,
  data_list$permuted = rbind(net_NA, net_notNA)

  return(data_list)
}

observedVSpermuted = function(data_list){
  data_list$result = data_list$permuted[data_list$observed, on = nodes$nodeX, allow.cartesian=TRUE]
  data_list$result[, mean(observed_statistic <= permuted_statistic), by = .(eval(nodes_call$nodeX), eval(nodes_call$nodeZ))]




  splitX = split(data_list$permuted, data_list$permuted$IDs_interactor_viral)
  splitX$U5TQE9

  pval_list = lapply(splitX, function(one_fold_enrichment_dist, observed){
    merged = one_fold_enrichment_dist[observed, nomatch = 0, on = "IDs_interactor_viral", allow.cartesian = T]
    #merged[, Pval := mean(observed_statistic <= permuted_statistic), by = IDs_domain_human]
    #unique(merged[,.(IDs_interactor_viral, IDs_domain_human, fold_enrichment, Pval)])
  }, data_list$observed)
  pval_list$U5TQE9
  pval_table = Reduce(rbind, pval_list)
  setorder(X, IDs_interactor_viral, permuted_statistic,IDs_interactor_human,IDs_interactor_viral_degree,IDs_domain_human,domain_count,observed_statistic)
  setorder(pval_table, IDs_interactor_viral, permuted_statistic,IDs_interactor_human,IDs_interactor_viral_degree,IDs_domain_human,domain_count,observed_statistic)
  all.equal(X, pval_table)
}
