##' download PSICQUIC query search result (only IMEx databases) to R library - queryPSICQUIC package - data directory
##' @name queryPSICQUICrlib
##' @author Vitalii Kleshchevnikov
##' @description saves PSIQUIC search results to the local R library or to the directory specified
##' @details \code{queryPSICQUICrlib} queries IMEx databases using PSICQUIC if the local copy of the query result is not available or if there is a new IntAct release. The functions will still work for non-IMEx databases but it will check the IntAct release date and that is not meaningful for tracking updates of other resources.
##' @inheritDotParams queryPSICQUIC -file
##' @param directory character, directory (non-existing directory will be created) where to store the result. The default is to store the result in R library inside the MItools package
##' @param releaseORdate character, if data has already been downloaded: which IntAct release or download date to read
##' @param just_list_releases logical, just list IntAct releases or download dates that are locally available for this query (format: try 2017Jul12 (IntAct release) or 20170712 (download date))
##' @return object of class "RAW_MItab25" or "RAW_MItab27" (list) containing 2 elements: data which is the query result (data.table) and metadata, a data.table which contains the query, the filepath, the format, the date or the IntAct release date, and the number of interactions per database
##' @import data.table
##' @importFrom R.utils gzip
##' @importFrom R.utils gunzip
##' @export queryPSICQUICrlib
##' @export print.RAW_MItab25
##' @export print.RAW_MItab27
##' @examples
##' {
##' queryPSICQUICrlib(query = "id:P74565 AND detmethod:\"MI:0018\"", format = "tab27", database = "imex")
##' }

queryPSICQUICrlib = function(..., directory = NULL, releaseORdate = NULL, just_list_releases = F){
  # to just list available releases for this query
  if(just_list_releases){
    query = paste0(...)
    filename = paste0("query_", query, ".tsv.gz")
    filename = gsub(":",".", filename)
    filename = gsub("\"",".", filename)
    filename = gsub(" ",".", filename)
    if(is.null(directory)){
      pkg_dir = paste0(.libPaths(), "/MItools", "/data/")[1]
      release_list = listReleases(pkg_dir, filename, query)
      print(release_list)
    } else {
      release_list = listReleases(directory, filename, query)
      print(release_list)
    }
    return(release_list)
  } else { # do the query and/or read the result
    # create a directory path name that contains information about either the date or the IntActRelease
    args = list(...)
    database = args$database
    if(is.null(directory)){
      pkg_dir = paste0(.libPaths(), "/MItools", "/data/")[1]
      # create data directory in /default.library/queryPSICQUIC/ if it doesn't exist
      if(!dir.exists(pkg_dir)) dir.create(pkg_dir)
      # find out last release date if the database is IntAct or imex, use query date if else, generate DirName
      pkg_dir_last_release = generateDirName(database, releaseORdate, pkg_dir)
    } else {
      pkg_dir_last_release = generateDirName(database, releaseORdate, directory)
    }

    # create directory for the last release date
    if(is.null(releaseORdate)) {
      if(!dir.exists(pkg_dir_last_release)) dir.create(pkg_dir_last_release)
    } else {
      if(!dir.exists(pkg_dir_last_release)) stop(paste0("no data for IntAct release or date: ", releaseORdate," in the directory: ", directory))
    }
    stopifnot(dir.exists(pkg_dir_last_release))
    # create a filepath name that uniquely identifies the search result
    filepath = paste0(pkg_dir_last_release, "query_", ..., ".tsv")
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

    if(is.null(releaseORdate)){
      suppressMessages({res[, date_time := lastIntActRelease()]})
    } else {
      res[, date_time := releaseORdate]
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

generateDirName = function(database, releaseORdate, directory){
  # find out last release date if the database is IntAct or imex, use query date if else
  if(is.null(database) | mean(database %in% "imex") == 1 | mean(database %in% "IntAct") == 1 | mean(database %in% "IntActFTP") == 1){
    if(is.null(releaseORdate)) {
      pkg_dir_last_release = paste0(directory, "IntActRelease_",lastIntActRelease(), "/")
    } else {
      pkg_dir_last_release = paste0(directory, "IntActRelease_", releaseORdate, "/")
    }
  } else {
    if(is.null(releaseORdate)) {
      date = Sys.Date()
      date = gsub("-","",date)
      pkg_dir_last_release = paste0(directory, "DownloadDate_",date, "/")
    } else {
      pkg_dir_last_release = paste0(directory, "DownloadDate_", releaseORdate, "/")
    }
  }
  return(pkg_dir_last_release)
}

listReleases = function(directory, filename, query){
  if(!dir.exists(directory)) stop(paste0("directory: ", directory, " doesn't exist"))
  release_list = list.files(directory, pattern = filename, recursive = T)
  if(length(release_list) == 0) {
    stop(paste0("no data for the query (",query,") is available in the directory:", directory))
  } else {
    release_list = gsub(paste0("/",filename), "", release_list)
    release_list = strsplit(release_list, "_")
  }
  return(release_list)
}
