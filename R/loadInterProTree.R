##' Load the hierarchy of relationships between InterPro's entries
##' @name loadInterProTree
##' @author Vitalii Kleshchevnikov
##' @param filename character, where to download and / or where to find InterPro Entry relationships tree file. Details: \link{https://www.ebi.ac.uk/interpro/download.html}
##' @param tree object of class "InterProTree", the output of \code{loadInterProTree} function
##' @param ParentChildTreeLines argument for non-exported \code{loadChildrenForLevel}: the lines of the InterPro Entry relationships tree file (character)
##' @param levels argument for non-exported \code{loadChildrenForLevel}: the mapping between the levels and the lines in ParentChildTreeLines
##' @param level argument for non-exported \code{loadChildrenForLevel}: integer, which level to get children for
##' @param regex argument for non-exported \code{loadChildrenForLevel}: regex that matches InterProID (IPR[[:digit:]]{6})
##' @param allchildren argument for non-exported \code{loadChildrenForLevel}: whether to return all children or only a specific level
##' @return object of S3 class "InterProTree" containing data.table that provides the mapping between layers, the lines of the InterPro Entry relationships tree file (character) and the mapping between the levels and the lines in ParentChildTreeLines
##' @seealso \link[MItools]{readInterProGFF3}
##' @export loadInterProTree
##' @export getLevelXchildren
##' @export print.InterProTree
##' @examples
##' InterProScan_domains_nonred = collapseByInterProID(InterProScan_features = InterProScan_domains, id_col = "Dbxref")
loadInterProTree = function(filename) {
  if(!file.exists(filename)) download("ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt", filename)
  ParentChildTreeLines = readLines(filename)

  levels = list()
  IPRpattern = paste0("IPR[[:digit:]]{6}")
  count = 1
  while (length(grep(paste0("^", IPRpattern), ParentChildTreeLines)) > 0) {
    levels[count][[1]] = grep(paste0("^", IPRpattern), ParentChildTreeLines)
    IPRpattern = paste0("--",IPRpattern)
    count = count + 1
  }

  level_tables = list()
  for (level in 1:(length(levels)-1)) {
    level_tables[[paste0("level", level)]] = loadChildrenForLevel(ParentChildTreeLines, levels, level = level, allchildren = F)
  }

  level_table = level_tables$level1
  for (level in 2:(length(levels)-1)) {
    lev_char = paste0("level",level)
    level_table = merge(level_table, level_tables[lev_char][[1]], by = lev_char, allow.cartesian = T, all.y = F, all.x = T)
  }
  col.order = colnames(level_table)[order(colnames(level_table))]
  level_table = unique(level_table[, col.order, with = F])
  level_table[!is.na(level5),]

  tree = list(level_table = level_table, levels2ParentChildTreeLines = levels, ParentChildTreeLines = ParentChildTreeLines)
  class(tree) = "InterProTree"
  return(tree)
}

getLevelXchildren = function(tree, level){
  loadChildrenForLevel(tree$ParentChildTreeLines, tree$levels2ParentChildTreeLines, level, allchildren = T)
}

loadChildrenForLevel = function(ParentChildTreeLines, levels, level = 1, regex = "IPR[[:digit:]]{6}", allchildren = F) {
  if(level == length(levels)) stop(paste0("This is the last level: ", level))
  tree = list()
  for (i in 1:length(levels[[level]])) {
    domainID = ParentChildTreeLines[levels[[level]][i]]
    str_start = gregexpr(regex, domainID)[[1]][1]
    level1name = substr(domainID, str_start, str_start + 8)
    count = levels[[level]][i] + 1
    domainIDs = character()
    if(i < length(levels[[level]])){
      inner = paste0("levels[[",1:level,"]]", collapse = ",")
      eval(parse(text = paste0("parent_level = c(", inner,")")))
      while (!count %in% parent_level) {
        domainID2 = ParentChildTreeLines[count]
        str_start = gregexpr(regex, domainID2)[[1]][1]
        domainIDs = c(domainIDs, substr(domainID2, str_start, str_start + 8))
        count = count + 1
      }
    }
    tree[[level1name]] = domainIDs
  }
  tree = data.table(parent = rep(names(tree), times = sapply(tree, length)), children = unlist(tree))
  setnames(tree, colnames(tree), c(paste0("level", level), "allchildren"))
  if(!allchildren){
    levelplus1 = ParentChildTreeLines[levels[[level + 1]]]
    str_start = unlist(gregexpr(regex, levelplus1))
    levelplus1 = substr(levelplus1, str_start, str_start + 8)
    tree = tree[allchildren %in% levelplus1,]
    setnames(tree, colnames(tree), c(paste0("level", level), paste0("level", level + 1)))
  }
  return(tree)
}

print.InterProTree = function(InterProTree) {
  cat("\n\n this object describes the hierarchy of relationships between InterPro's entries (i.e. families and their subfamilies) \n\n")
  cat("details: https://www.ebi.ac.uk/interpro/download.html \n\n")
  cat("\n these relationships are stored in $level_table \n")
  print(InterProTree$level_table)
  cat("\n if you want all children of for specific level use getLevelXchildren function on this object: \n")
  cat("\n getLevelXchildren(tree = this_object, level = 1) \n")
}
