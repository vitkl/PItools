##' functions that are called inside permutationPval
##' @name permutationPvalHelper
##' @rdname permutationPvalHelper
##' @export permutationPvalHelper
permutationPvalHelper = NULL
##' extract and match attribute column names to data.tables that contain nodes which those attributes describe
##' @name formula2colnames
##' @rdname permutationPvalHelper
##' @param inheritParams permutationPval
##  node_attr interactions2permute associations2test select_nodes also_permuteYZ
##' @param cols list of column names
##' @param nodes list of node names (column names for columns that contain node names)
##' @param nodes_call list of node names, each class call (to be evaluated within data.table expressions DT[,eval(nodes_call$nodeX)])
##' @param data_list list that contains data.tables necessary for these functions to work
##' @param by_cols modified left-hand side of \code{statistic} (\code{\link[MItools]{permutationPval}}), class call
##' @param exprs right-hand side of \code{statistic} (\code{\link[MItools]{permutationPval}}), class call, \code{exprs} is evaluated by \code{by_cols}; in data.table synthax: DT[, statistic := eval(exprs), by = .(eval(by_cols))]
##' @param includeAssociations logical, if calculating statistic requires columns that contain attribute of X-Z pair associations table will be merged
##' @details \code{formula2colnames} extracts nodes, extracts attributes, attaches attributes to a table if all nodes are present in that table and in case of associations table both varibles in a formula should be nodeX and nodeZ
##' @return \code{formula2colnames} returns list of column names for XY, YZ and association tables
##' @import data.table
##' @import BiocGenerics
##' @author Vitalii Kleshchevnikov
##' @usage cols = formula2colnames(node_attr, cols, nodes)
formula2colnames = function(node_attr, cols, nodes){
  # if not formula is provided
  if(!is.formula(node_attr)) stop(paste("node_attr is provided but is not a formula: class - ", class(node_attr), "; content - ",paste(paste0("[[",1:length(node_attr),"]]"), node_attr, collapse = " ")))

  # extract nodes
  node_temp = all.vars(node_attr[[2]])
  if(mean(node_temp %in% unlist(nodes)) != 1) stop(paste0("node_attr is supplied for node: ",node_temp," , however, one or more of these nodes are not supplied in interactions2permute or associations2test"))
  # extract attributes
  attribute_temp = all.vars(node_attr[[3]])
  # attach attributes to a table if all nodes are present in that table
  if(mean(node_temp %in% c(nodes$nodeX, nodes$nodeY)) == 1) cols$interactionsXY_cols = c(cols$interactionsXY_cols, attribute_temp)
  if(mean(node_temp %in% c(nodes$nodeY, nodes$nodeZ)) == 1) cols$interactionsYZ_cols = c(cols$interactionsYZ_cols, attribute_temp)
  # and in case of associations table both varibles in a formula should be nodeX and nodeZ
  if(mean(c(nodes$nodeX, nodes$nodeZ) %in% node_temp) == 1) cols$associations_cols = c(cols$associations_cols, attribute_temp)
  return(cols)
}

##' @name node_attr2colnames
##' @rdname permutationPvalHelper
##' @return \code{node_attr2colnames} returns the same as \code{formula2colnames} but when node_attr2colnames is a list of formulas
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
##' @return \code{formula2nodes} returns list of 3 node column names extracted from formulas supplied through \code{interactions2permute}, \code{associations2test}
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

##' @name checkSelectNodes
##' @rdname permutationPvalHelper
##' @details \code{checkSelectNodes} checks if you want to select nodes based on the attributes of those nodes + splits a formula which selects a combination of nodes based on their shared properties (like count of X and Z pairs) into a list of formulas each selecting separate nodes
##' @return \code{checkSelectNodes} returns a formula or a list of formulas if they have passed the check - error if otherwise
##' @usage select_nodes = checkSelectNodes(select_nodes, node_attr, nodes, cols)
checkSelectNodes = function(select_nodes, node_attr, nodes, cols){

  if(!(is.formula(select_nodes))) stop("select_nodes is not formula")

  select_node = all.vars(select_nodes[[2]])
  select_attribute = all.vars(select_nodes[[3]])
  all_nodes = unlist(nodes)
  all_cols = unique(unlist(cols))

  if(length(select_node) >= 3) stop("you cannot select 3 or more nodes based on their shared attribute, try filtering the data before using permutationPval")
  if(mean(select_node %in% all_nodes) != 1) stop(paste0("select_nodes is supplied for node(s): ",paste0(select_node, collapse = " and ")," , however, one or more of these nodes are not supplied in interactions2permute or associations2test"))
  if(mean(select_attribute %in% all_cols) != 1) stop(paste0("select_nodes is supplied with attribute(s): ",paste0(select_attribute, collapse = " and ")," , however, one or more of these attributes are not supplied in node_attr"))

  if(is.null(node_attr)){
    if(!((mean(select_node %in% all_nodes) == 1) &
       (mean(select_attribute %in% all_nodes) == 1))) stop(paste0("select_nodes is supplied for node(s): ",paste0(select_node, collapse = " and ")," , however, one or more of these nodes are not supplied in interactions2permute or associations2test"))
  } else {
    if(!(is.formula(node_attr))) stop("node_attr is not formula")
    node = all.vars(node_attr[[2]])
    node_attribute = all.vars(node_attr[[3]])

    # if select_nodes select nodes by attribute in node_attr
    if((mean(select_node %in% node) == 1) & (mean(node %in% select_node) == 1)){
      if(mean(select_attribute %in% c(all_nodes, node_attribute)) != 1) stop(paste0("select_nodes asks to select ",paste0(select_node, collapse = " and ")," by attribute (",paste0(select_attribute, collapse = " "),") not specified in node_attr as the attribute of ",paste0(select_node, collapse = " and ")))
    }

  }

  if(length(select_node) == 1) return(select_nodes)
  if(length(select_node) == 2) {
    select_node_list = list()
    select_node_list[[1]] = as.formula(paste0(select_node[1],"~", as.expression(select_nodes[[3]])))
    select_node_list[[2]] = as.formula(paste0(select_node[2],"~", as.expression(select_nodes[[3]])))
    return(select_node_list)
  }
}

##' @name checkSelectNodesList
##' @rdname permutationPvalHelper
##' @details \code{checkSelectNodesList} calls \code{checkSelectNodes} for each element in a list (or formula) \code{select_nodes} and each element in a list (or formula) \code{node_attr}. If \code{select_nodes} or \code{node_attr} is a formula then it is converted to a list of 1L
##' @return \code{checkSelectNodesList} returns a list of formulas (even of lenght 1) if they have passed the check - error if otherwise
##' @usage select_nodes = checkSelectNodesList(select_nodes, node_attr, nodes, cols)
checkSelectNodesList = function(select_nodes, node_attr, nodes, cols){
  if(is.formula(select_nodes)) select_nodes = list(select_nodes)
  select_nodes_new = list()
  if(is.formula(node_attr)) node_attr = list(node_attr)
  for (select_formula in select_nodes) {
    for (node_attr_formula in node_attr) {
      select_nodes_temp = checkSelectNodes(select_formula, node_attr_formula, nodes, cols)
      select_nodes_new = c(select_nodes_new, select_nodes_temp)
    }
  }
  select_nodes = unique(select_nodes_new)
}

##' @name filterByFormula
##' @rdname permutationPvalHelper
##' @details \code{filterByFormula} filters one or 2 of the XY, YZ or association tables based on a condition defining which nodes to keep (only one node type: X, Y or Z)
##' @return \code{filterByFormula} returns data_list in which interactionsXY, interactionsYZ, associations data.table-s were filtered by formula
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

##' @name filterByList
##' @rdname permutationPvalHelper
##' @details \code{filterByList} filters one or 2 of the XY, YZ or association tables based on a condition defining which nodes to keep (only one node type: X, Y or Z), takes a list of formulas or a formula as an argument
##' @return \code{filterByList} returns data_list in which interactionsXY, interactionsYZ, associations data.table-s were filtered by each formula if select_nodes was a list
##' @usage data_list = filterByList(data_list, select_nodes, cols, nodes)
filterByList = function(data_list, select_nodes, cols, nodes){
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
  } else if(is.null(select_nodes)) NULL else stop("select_nodes is provided but is neither a list nor a formula")
  return(data_list)
}

##' @name calcObservedStatistic
##' @rdname permutationPvalHelper
##' @details \code{calcObservedStatistic} calculates Observed Statistic by evaluating \code{exprs} by \code{by_cols}, also calculates how many nodes Y are missing node Z for each node X, saves result to \code{data_list$observed}. includeAssociations specifies if \code{exprs} required any attributes of both X and Z.
##' @return \code{calcObservedStatistic} returns \code{data_list}
##' @usage data_list = calcObservedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations)
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

##' @name calcPermutedStatistic
##' @rdname permutationPvalHelper
##' @details \code{calcPermutedStatistic} calculates Permuted Statistic by permuting node Y column in XY (and, optionally in YZ, if \code{also_permuteYZ} is true) and calculating statistic by evaluating \code{exprs} by \code{by_cols}, unlike \code{calcObservedStatistic} doesn't calculate how many nodes Y are missing node Z for each node X. Saves result to \code{data_list$permuted}.
##' @return \code{calcPermutedStatistic} returns \code{data_list}
##' @usage data_list = calcPermutedStatistic(data_list, by_cols, exprs, nodes, nodes_call, includeAssociations, also_permuteYZ)
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

##' @name observedVSpermuted
##' @rdname permutationPvalHelper
##' @details \code{observedVSpermuted} compares Observed Statistic to Permuted Statistic, records how many times permuted statistic is at least as large as the observed statistic (\code{higher_counts}) and how many of X-Z pairs have non-NA Observed Statistic or Permuted Statistic (\code{not_missing}. In each permutation round, node X might have many nodes Z, so when calculating empirical p-value all these comparisons should be counted (not just the number of permutations). On the other hand, not in all permutation rounda specific observed XZ pair might occur, so we need to cound non-NA "Observed Statistic to Permuted Statistic" comparisons
##' @return \code{observedVSpermuted} returns data.table contatining the following columns: "higher_counts", "not_missing", nodeX names, nodeZ names
##' @usage data_list_temp = MItools:::observedVSpermuted(data_list, nodes_call, nodes)
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

##' @name aggregatePermutations
##' @rdname permutationPvalHelper
##' @details \code{aggregatePermutations} aggregates the result of \code{observedVSpermuted} from many permutation rounds (output of replicate, class matrix)
##' @return \code{aggregatePermutations} returns the data.table of the same format as \code{observedVSpermuted}
##' @usage temp2_inner = MItools:::aggregatePermutations(temp_inner, nodes, nodes_call)
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
