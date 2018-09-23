##' parReplicate function from: https://stackoverflow.com/questions/19281010/simplest-way-to-do-parallel-replicate
##' @details arguments as in \code{\link{parSapply}} and \code{\link{replicate}}
##' @importFrom BiocGenerics parSapply
##' @export parReplicate
parReplicate = function(cl, n, expr, simplify=TRUE, USE.NAMES=TRUE){
  parSapply(cl, integer(n), function(i, ex) eval(ex, envir=.GlobalEnv),
            substitute(expr), simplify=simplify, USE.NAMES=USE.NAMES)
}
