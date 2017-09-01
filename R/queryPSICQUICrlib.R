##' download PSICQUIC query search result (only IMEx databases) to R library - queryPSICQUIC package - data directory
##' @name queryPSICQUICrlib
##' @author Vitalii Kleshchevnikov
##' @description saves PSIQUIC search results to the local R library or to the directory specified
##' @details \code{queryPSICQUICrlib} queries IMEx databases using PSICQUIC if the local copy of the query result is not available or if there is a new IntAct release. The functions will still work for non-IMEx databases but it will check the IntAct release date and that is not meaningful for tracking updates of other resources.
##' @inheritDotParams queryPSICQUIC -file
##' @param directory character, directory (non-existing directory will be created) where to store the result. The default is to store the result in R library inside the MItools package
##' @return object of class "RAW_MItab25" or "RAW_MItab27" (list) containing 2 elements: data which is the query result (data.table) and metadata, a data.table which contains the query, the filepath, the format, the date or the IntAct release date, and the number of interactions per database
##' @import data.table
##' @importFrom R.utils gzip
##' @importFrom R.utils gunzip
##' @export queryPSICQUICrlib
##' @export print.RAW_MItab25
##' @export print.RAW_MItab27
##' @examples MItools:::downloadPSICQUICrlib(query = "id:P74565 AND detmethod:\"MI:0018\"", format = "tab27", database = "imex")

queryPSICQUICrlib = function(..., directory = NULL){
  # create data directory in /default.library/queryPSICQUIC/ if it doesn't exist
  if(is.null(directory)){
    pkg_dir = paste0(.libPaths(), "/MItools", "/data/")[1]
    if(!dir.exists(pkg_dir)) dir.create(pkg_dir)
    # create a filepath that uniquely identifies the search result and contains a timestamp
    # find out last release date if the database is IntAct or imex, use query date if else
    args = list(...)
    database = args$database
    if(is.null(database) | mean(database %in% "imex") == 1 | mean(database %in% "IntAct") == 1){
      pkg_dir_last_release = paste0(pkg_dir, lastIntActRelease())
    } else {
      pkg_dir_last_release = paste0(pkg_dir, Sys.Date())
    }
  } else {
    # find out last release date if the database is IntAct or imex, use query date if else
    args = list(...)
    database = args$database
    if(is.null(database) | mean(database %in% "imex") == 1 | mean(database %in% "IntAct") == 1){
      pkg_dir_last_release = paste0(directory, lastIntActRelease())
    } else {
      pkg_dir_last_release = paste0(directory, Sys.Date())
    }
    }

  # create directory for the last release date
  if(!dir.exists(pkg_dir_last_release)) dir.create(pkg_dir_last_release)
  stopifnot(dir.exists(pkg_dir_last_release))
  filepath = paste0(pkg_dir_last_release, "/query_", ..., ".tsv")
  filepath = gsub(":",".", filepath)
  filepath = gsub("\"",".", filepath)
  filepath = gsub(" ",".", filepath)
  filepath.gz = paste0(filepath, ".gz")

  # read results if the data from the latest release has been downloaded
  if(file.exists(filepath.gz)){
    message("found local copy of the data from the latest release ... reading into R")
    res = fread(paste0(filepath,"log"), header = T, stringsAsFactors = F, sep = "\t")
  }
  # download results using queryPSICQUIC function if the data from the latest release hasn't been downloaded
  if(!file.exists(filepath.gz)){
    message("downloading using PSICQUIC")
    res = queryPSICQUIC(file = filepath, ...)
    fwrite(res, file = paste0(filepath,"log"), sep = "\t")
    R.utils::gzip(filename = filepath,
                  destname = filepath.gz,
                  overwrite = T, skip = F, remove = T)
  }

  if(is.null(directory)){
    res[, date_time := gsub(pkg_dir, "",pkg_dir_last_release)]
  } else {
    res[, date_time := gsub(directory, "",pkg_dir_last_release)]
    }
  # read the result
  R.utils::gunzip(filename = filepath.gz,
                  destname = filepath,
                  overwrite = T, skip = F, remove = F)
  result = list(data = fread(filepath, header = T, stringsAsFactors = F), metadata = res)
  unlink(filepath)
  # add class attr.: 1. RAW means mitab as is, 2. format: MITAB25 or MITAB27
  class(result) = paste0("RAW_MI",unique(res$format)[1])

  return(result)
}

#print methods
print.RAW_MItab25 = function(data){
  cat(paste0("\n` Object of class RAW_MItab25, for query, file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
print.RAW_MItab27 = function(data){
  cat(paste0("\n` Object of class RAW_MItab27, for query, file, format, databases, date: `\n"))
  print(data$metadata)
  cat("\n` view of the $data: `\n")
  print(data$data)
}
