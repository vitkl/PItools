##' download PSICQUIC query search result (only IMEx databases) to R library - queryPSICQUIC package - data directory
##' @name queryPSICQUICrlib
##' @author Vitalii Kleshchevnikov
##' @description saves PSIQUIC search results to the local R library or to the directory specified
##' @details \code{queryPSICQUICrlib} queries IMEx databases using PSICQUIC if the local copy of the query result is not available or if there is a new IntAct release. The functions will still work for non-IMEx databases but it will check the IntAct release date and that is not meaningful for tracking updates of other resources.
##' @param ... args passed to \code{\link{queryPSICQUIC}}
##' @return data.table containing the query result
##' @import data.table
##' @examples queryPSICQUIC:::downloadPSICQUICrlib(query = "id:P74565 AND detmethod:\"MI:0018\"", format = "tab27", database = "imex")

queryPSICQUICrlib = function(..., directory = NULL){
  # create data directory in /default.library/queryPSICQUIC/ if it doesn't exist
  if(is.null(directory)){
    pkg_dir = paste0(.libPaths(), "/MItools", "/data/")[1]
    if(!dir.exists(pkg_dir)) dir.create(pkg_dir)
    # create a filepath that uniquely identifies the search result and contains a timestamp
    # find out last release date
    pkg_dir_last_release = paste0(pkg_dir, lastIntActRelease())
  }
  if(!is.null(directory)) pkg_dir_last_release = paste0(directory, lastIntActRelease())

  # create directory for the last release date
  if(!dir.exists(pkg_dir_last_release)) dir.create(pkg_dir_last_release)
  stopifnot(dir.exists(pkg_dir_last_release))
  filepath = paste0(pkg_dir_last_release, "/query_", ..., ".tsv")
  filepath = gsub(":",".", filepath)
  filepath = gsub("\"",".", filepath)
  filepath = gsub(" ",".", filepath)

  # download results using queryPSICQUIC function if the data from the latest release haven't been downloaded
  if(file.exists(filepath)){
    message("found local copy of the data from the latest release ... reading into R")
    res = fread(paste0(filepath,"log"), header = T, stringsAsFactors = F, sep = "\t")
    print(res)
  }
  if(!file.exists(filepath)){
    message("downloading using PSICQUIC")
    res = queryPSICQUIC(file = filepath, ...)
    fwrite(res, file = paste0(filepath,"log"), sep = "\t")
    print(res)
  }

  # read the result and return
  return(fread(filepath, header = T, stringsAsFactors = F))

}


