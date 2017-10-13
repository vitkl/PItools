##' Generate a QSLIMFinder command
##' @rdname QSLIMFinderCommand
##' @name QSLIMFinderCommand
##' @author Vitalii Kleshchevnikov
##' @param file_list data.table containing path to files and directories for QSLIMFinder: fastafile, queryfile, outputdir, outputfile
##' @param i integer, which set of files and directories to choose from \code{file_list}
##' @param slimpath relative path (from the project folder) to qslimfinder.py
##' @param blast relative path (from the project folder) to /bin/ folder in BLAST package
##' @param iupred relative path (from the project folder) to iupred (compiled executable)
##' @param options any options from QSLIMFinder
##' @param LSF_cluster_mode logical, if FALSE \code{LSF_cluster_par} and \code{LSF_project_path} are ignored
##' @param LSF_cluster_par a string that will launch LSF job
##' @param LSF_project_path absolute path on LSF cluster
##' @return character vector (1L), bash command that will lauch QSLIMFinder locally or as a job on LSF cluster
##' @import data.table
##' @export QSLIMFinderCommand
##' @seealso \code{\link{listInteractionSubsetFASTA}}, \code{\link{runQSLIMFinder}}
QSLIMFinderCommand = function(file_list, i = 1,
                              slimpath = "../software/slimsuite/tools/qslimfinder.py",
                              blast = "../software/ncbi_blast_2.6.0/bin/",
                              iupred = "../software/iupred/iupred",
                              options = "dismask=T consmask=F cloudfix=F probcut=0.1",
                              LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                              LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/")
{
  dirs = file_list[i,]

  dirs_command = paste0("resdir=", LSF_project_path, dirs$outputdir,
                        " resfile=", LSF_project_path, dirs$outputfile,
                        " seqin=", LSF_project_path, dirs$fastafile,
                        " query=", LSF_project_path, dirs$queryfile," ")
  slimpath = paste0(LSF_project_path,slimpath)
  blast = paste0("blast+path=",LSF_project_path, blast)
  iupred = paste0("iupath=",LSF_project_path, iupred)

  command = paste(LSF_cluster_par, "python", slimpath, blast, iupred, options, dirs_command)
  #system(command, ignore.stdout = T, ignore.stderr = T)
  return(command)
}

##' @rdname QSLIMFinderCommand
##' @name mQSLIMFinderCommand
##' @import data.table
##' @export mQSLIMFinderCommand
##' @param log_dir character, directory where to write log files (stout and sterr)
##' @param write_log whether to save stout and sterr
##' @return list containing: 1. command to set up enviromental variable IUPred_PATH; 2. character vector of bash commands that will lauch QSLIMFinder as a job on LSF cluster; 3, 4, 5 - directories where LSF should write stout and sterr
##' @examples
##' all_commands = mQSLIMFinderCommand(file_list = forSLIMFinder_file_list,
##'     slimpath = "../software/cluster/slimsuite/tools/qslimfinder.py",
##'     blast = "../software/cluster/ncbi_blast_2.6.0/bin/",
##'     iupred = "../software/cluster/iupred/iupred",
##'     options = "dismask=T consmask=F cloudfix=T probcut=0.3 iuchdir=T",
##'     LSF_cluster_mode = T,
##'     LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
##'     LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"")
mQSLIMFinderCommand = function(file_list,
                               slimpath = "../software/slimsuite/tools/qslimfinder.py",
                               blast = "../software/ncbi_blast_2.6.0/bin/",
                               iupred = "../software/iupred/iupred",
                               options = "dismask=T consmask=F cloudfix=F probcut=0.5",
                               LSF_cluster_mode = F,
                               LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                               LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 200 -R \"rusage[mem=200]\"",
                               log_dir = "./SLIMFinder/log_dir/",
                               write_log = T, recursive = F)
{
  if(LSF_cluster_mode){
    commands1 = paste0("export IUPred_PATH=",LSF_project_path,iupred)

    if(write_log){
      log_dirfull = paste0(LSF_project_path, log_dir)
      log_dirlog = paste0(LSF_project_path, log_dir, "log/")
      log_direrror = paste0(LSF_project_path, log_dir, "error/")
    }

    commands = sapply(1:nrow(file_list), function(i){
      name = paste0("interactors_of",file_list$interactors_of[i],"QSLIMFinder_query", file_list$QSLIMFinder_query[i])
      if(write_log){
        LSF_cluster_par = paste0(LSF_cluster_par," -o ", log_dirlog, name,"log -e ", log_direrror, name, "error")
      } else {
        LSF_cluster_par = paste0(LSF_cluster_par," -o /dev/null -e /dev/null")
      }

      QSLIMFinderCommand(file_list = file_list[i,], slimpath = slimpath,
                         blast = blast, iupred = iupred, options = options,
                         LSF_cluster_par = LSF_cluster_par,
                         LSF_project_path = LSF_project_path)
    })

  } else {
    commands1 = paste0("export IUPred_PATH=",iupred)
    commands = sapply(1:nrow(file_list), function(i){
      QSLIMFinderCommand(file_list = dirs, slimpath = slimpath,
                         blast = blast, iupred = iupred, options = options,
                         LSF_cluster_par = "",
                         LSF_project_path = "")
    })
  }

  if(write_log){
    return(list(set_env_var = commands1, run = commands,
                log_dirfull = log_dirfull, log_dirlog = log_dirlog, log_direrror = log_direrror))
  } else {
    return(list(set_env_var = commands1, run = commands,
                log_dirfull = NULL, log_dirlog = NULL, log_direrror = NULL))
  }
}
