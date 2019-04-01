##' load all data from IntAct ftp
##' @name loadIntActFTP
##' @author Vitalii Kleshchevnikov
##' @description \code{loadIntActFTP} loads all data stored on \link{ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/} in intact.txt file. Details: \link{https://www.ebi.ac.uk/intact/downloads}. This file contains the data from the following databases "IntAct", "MINT", "DIP", "bhf-ucl", "MPIDB", "MatrixDB", "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo"
##' @param dir directory where to save/look for the local copy
##' @param release which locally saved IntAct release to load (the default is to load the latest and read it into R)
##' @return \code{loadIntActFTP} saves intact.txt to a file, returns the object of class RAW_MItab27
##' @import data.table
##' @importFrom R.utils gzip
##' @importFrom R.utils gunzip
##' @importFrom downloader download
##' @export loadIntActFTP
##' @examples
##' {
##' loadIntActFTP("./", release = NULL)
##' }
##' @author Vitalii Kleshchevnikov
loadIntActFTP = function(dir, release = NULL){
  if(is.null(release)) file_name = paste0(dir,"intact",lastIntActRelease(), ".txt") else file_name = paste0(dir,"intact",release, ".txt")
  file_name.gz = paste0(file_name, ".gz")
  file_name.zip = gsub("txt", "zip", file_name)

  # download MI-TAB 2.7 from IntAct ftp
  if(!file.exists(file_name.gz)) {
    if(!is.null(release)) stop(paste0("no data for IntAct release", release," in the directory: ", dir, ", set release to NULL do load the latest release"))
    message("... dowloading from IntAct ftp ...")
    download(url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip", destfile = file_name.zip)
    unzip(file_name.zip, exdir = dir)
    unlink(paste0(dir, "intact_negative.txt"))
    unlink(file_name.zip)
    gzip(filename = paste0(dir, "intact.txt"), destname = file_name.gz,
         remove = T, overwrite = T)

  } else {message("... loading local copy ...")}

  gunzip(filename = file_name.gz, destname = file_name, remove = F, overwrite = T)
  intact = fread(file_name, header = T, stringsAsFactors = F)
  unlink(file_name)

  result = list(data = intact, metadata = "This object contains the data from ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt and contains all molecular interaction data from the following databases: \"IntAct\", \"MINT\", \"DIP\", \"bhf-ucl\", \"MPIDB\", \"MatrixDB\", \"HPIDb\",\"I2D-IMEx\",\"InnateDB-IMEx\", \"MolCon\", \"UniProt\", \"MBInfo\"")

  class(result) = "RAW_MItab27"
  return(result)
}
