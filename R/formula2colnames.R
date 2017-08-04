##' extract and match attribute column names to data.tables that contain nodes which those attributes describe
##' @name formula2colnames
##' @aliases node_attr2colnames
##' @aliases formula2nodes
##' @aliases data_list
##' @param node_attr formula (formula2colnames and node_attr2colnames) or list of formulas (node_attr2colnames)
##' @param cols list of column names
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
formula2colnames = function(node_attr, cols, nodes){
  # if not formula is provided
  if(!is.formula(node_attr)) stop(paste("node_attr is provided but is not a formula: ", paste(class(node_attr),node_attr, collapse = " ")))

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
  # extract node column as expression
  node_ftemp = as.expression(select_nodes[[2]])
  # extract condition as expression
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
