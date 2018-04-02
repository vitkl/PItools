##' Run CompariMotif3
##' @rdname runCompariMotif3
##' @name runCompariMotif3
##' @author Vitalii Kleshchevnikov
##' @param input_file a path to a tsv file containing Name and Pattern columns
##' @param slimpath path to /slimsuite/tools/
##' @param dbpath path to directory where to save and keep ELM database (\link{http://elm.eu.org/}) or other database of linear motifs in a format required by comparimotif_V3: \link{http://rest.slimsuite.unsw.edu.au/docs&page=module:comparimotif_V3}
##' @param dburl url where to download database
##' @param parameters any commandline parameters for comparimotif_V3 as a single character string
##' @param out_file path to a file where comparimotif_V3 should save the result
##' @param LSF_project_path path to the project directory on the LSF cluster file system. If on LSF is to be used all other paths should be relative to this directory, alternatively, set this argument to empty character. Final path is \code{LSF_project_path} concatenated with \code{input_file}, \code{slimpath}, \code{out_file} or \code{dbpath}. When using this on a local machine set this argument to empty character.
##' @param run logical, run comparimotif_V3 (T) or just output the command (F), can be useful for check correcteness
##' @param with compare motifs in \code{input_file} to the database in \code{dbpath} ("db") or to \code{input_file} ("self")
##' @return system command to be run or that was run
##' @import data.table
##' @export runCompariMotif3
##' @seealso \code{\link{writeInteractionSubsetFASTA_list}}, \code{\link{runQSLIMFinder}}
runCompariMotif3 = function(input_file = "./SLIMFinder/result/motifs.txt",
                            slimpath = "../software/slimsuite/tools/",
                            dbpath = "./data_files/",
                            dburl = "http://elm.eu.org/elms/elms_index.tsv",
                            parameters = "unmatched=T motific=T",
                            out_file = "./SLIMFinder/result/comparimotif.tdt",
                            LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                            run = F, with = c("self", "db")[2]){
  if(with == "db"){
    elm_filename = paste0(LSF_project_path, dbpath, data.table::year(Sys.Date()), "elms_index.tsv")
    if(!file.exists(elm_filename)) download.file(dburl, elm_filename)
  }
  comparimotif_call = paste0("python ", LSF_project_path, slimpath, "comparimotif_V3.py")
  input_res = paste0("motifs=", LSF_project_path, input_file)
  if(with == "self"){
    input_database = ""
  } else if(with == "db"){
    input_database = paste0("searchdb=", elm_filename)
  } else stop("'with' is provided but is not 'self' or 'db'")
  out_file = paste0("resfile=", LSF_project_path, out_file)
  command = paste(comparimotif_call, input_res, input_database, parameters, out_file)
  if(run) system(command, wait = F, ignore.stdout = T)
  return(command)
}
