##' Generate a QSLIMFinder command
##' @rdname QSLIMFinderCommand
##' @name QSLIMFinderCommand
##' @author Vitalii Kleshchevnikov
##' @param file_list data.table containing path to files and directories for QSLIMFinder: fastafile, queryfile, outputdir, outputfile
##' @param i integer, which set of files and directories to choose from \code{file_list}
##' @param slimpath relative path (from the project folder) to the directory containing qslimfinder.py or slimfinder.py
##' @param blast relative path (from the project folder) to /bin/ folder in BLAST package
##' @param iupred relative path (from the project folder) to iupred (compiled executable)
##' @param options any options from QSLIMFinder
##' @param LSF_cluster_mode logical, if FALSE \code{LSF_cluster_par} and \code{LSF_project_path} are ignored
##' @param LSF_cluster_par a string that will launch LSF job
##' @param LSF_project_path absolute path on LSF cluster
##' @param analysis_type character 1L, qslimfinder or slimfinder. slimfinder doesn't need query: identical datasets with the same query will be removed and query files not written)
##' @return character vector (1L), bash command that will lauch QSLIMFinder locally or as a job on LSF cluster
##' @import data.table
##' @export QSLIMFinderCommand
##' @seealso \code{\link{listInteractionSubsetFASTA}}, \code{\link{runQSLIMFinder}}
QSLIMFinderCommand = function(file_list, i = 1,
                              slimpath = "../software/slimsuite/tools/",
                              blast = "../software/ncbi_blast_2.6.0/bin/",
                              iupred = "../software/iupred/iupred",
                              options = "dismask=T consmask=F cloudfix=F probcut=0.1",
                              LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                              LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                              analysis_type = c("qslimfinder", "slimfinder")[1])
{
  if(!analysis_type %in% c("qslimfinder", "slimfinder")) stop("analysis_type should be one of the \"qslimfinder\", \"slimfinder\"")
  dirs = file_list[i,]

  if(analysis_type == "qslimfinder"){
    dirs_command = paste0("resdir=", LSF_project_path, dirs$outputdir,
                          " resfile=", LSF_project_path, dirs$outputfile,
                          " seqin=", LSF_project_path, dirs$fastafile,
                          " query=", LSF_project_path, dirs$queryfile," ")
  }
  if(analysis_type == "slimfinder"){
    dirs_command = paste0("resdir=", LSF_project_path, dirs$outputdir,
                          " resfile=", LSF_project_path, dirs$outputfile,
                          " seqin=", LSF_project_path, dirs$fastafile, " ")
  }

  slimpath = paste0(LSF_project_path, slimpath, analysis_type, ".py")
  blast = paste0("blast+path=",LSF_project_path, blast)
  iupred = paste0("iupath=",LSF_project_path, iupred)

  command = paste(LSF_cluster_par, "python", slimpath, blast, iupred, options, dirs_command)
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
                               slimpath = "../software/slimsuite/tools/",
                               blast = "../software/ncbi_blast_2.6.0/bin/",
                               iupred = "../software/iupred/iupred",
                               options = "dismask=T consmask=F cloudfix=F probcut=0.5",
                               LSF_cluster_mode = F,
                               LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                               LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 200 -R \"rusage[mem=200]\"",
                               log_dir = "./SLIMFinder/log_dir/",
                               write_log = T, recursive = F,
                               analysis_type = c("qslimfinder", "slimfinder")[1])
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
                         LSF_project_path = LSF_project_path,
                         analysis_type = analysis_type)
    })

  } else {
    commands1 = paste0("export IUPred_PATH=",iupred)
    commands = sapply(1:nrow(file_list), function(i){
      QSLIMFinderCommand(file_list = dirs, slimpath = slimpath,
                         blast = blast, iupred = iupred, options = options,
                         LSF_cluster_par = "",
                         LSF_project_path = "",
                         analysis_type = analysis_type)
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

##' @rdname QSLIMFinderCommand
##' @name groupQSLIMFinderCommand
##' @import data.table
##' @export groupQSLIMFinderCommand
##' @param commands list returned by \code{mQSLIMFinderCommand()} containing: 1. command to set up enviromental variable IUPred_PATH; 2. character vector of bash commands that will lauch QSLIMFinder as a job on LSF cluster; 3, 4, 5 - directories where LSF should write stout and sterr
##' @param InteractionSubsetFASTA_list object of class InteractionSubsetFASTA_list containing: FASTA sequences for interacting proteins, molecular interaction data they correspond to. Each element of a list contains input for individual QSLIMFinder run.
##' @param sh_dir directory within dataset directory (dataset_name) where to write batch command .sh files
##' @param LSF_project_path full path to the project where dataset directory (dataset_name) is located.
##' @param dataset_name character, name of the dataset, such as "SLIMFinder" or "SLIMFinder_Vidal"
##' @param N_seq size of the batch. Groups QSLIMFinder jobs until the next job doesn't fit into N_seq, if jobs is larger than N_seq a single job will be written to a batch .sh file.
##' @param write_log FALSE will not allow runQSLIMFinder to detect crashed jobs
##' @return QSLIMFinderCommand split into batches by N_seq sequences. List containing: 1. command to set up enviromental variable IUPred_PATH; 2. character vector of bash commands that will lauch $SHELL as a job on LSF cluster and run QSLIMFinder commands from batch file; 3, 4, 5 - directories where LSF should write stout and sterr
groupQSLIMFinderCommand = function(commands, InteractionSubsetFASTA_list, sh_dir = "/sh_dir/", LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/", dataset_name = "SLIMFinder", N_seq = 200, write_log = T) {
  if(class(InteractionSubsetFASTA_list) != "InteractionSubsetFASTA_list") stop("InteractionSubsetFASTA_list should be of class \"InteractionSubsetFASTA_list\"")
  if(mean(c("set_env_var", "run", "log_dirfull", "log_dirlog", "log_direrror") %in% names(commands)) != 1) stop("commands should be the output of mQSLIMFinderCommand()")

  sh_dir = paste0("./", dataset_name,sh_dir)
  if(!dir.exists(sh_dir)) dir.create(sh_dir)

  seqNUM = sapply(InteractionSubsetFASTA_list$fasta_subset_list, function(x) length(unique(x)))
  cum_seqNUM = cumsum(seqNUM)
  batches = seq(N_seq, max(cum_seqNUM), N_seq)
  batches_ass = sapply(batches, function(batch, cum_seqNUM) {
    cum_seqNUM < batch
  }, cum_seqNUM)
  batches_ass = apply(batches_ass, 1, function(sample) which(sample)[1])
  batches_ass = paste0("batch_", batches_ass)
  seq_per_batch = sapply(split(seqNUM, batches_ass), sum)
  seq_per_batch = seq_per_batch[unique(names(seq_per_batch))]
  commands_run = split(commands$run, batches_ass)
  bsub = sapply(seq_along(commands_run), function(i) {
    bsub_temp = unique(gsub("python.+","",commands_run[[i]]))
    hps_dir = paste0(LSF_project_path, dataset_name,"/")

    batch_name = names(commands_run[i])

    sh_file = paste0(sh_dir, batch_name, ".sh")
    command_temp = unique(gsub("^.+ python","python",commands_run[[i]]))
    write(command_temp, sh_file)

    bsub_temp = unique(gsub("-o .+/log_dir/log/.+ -e .+/log_dir/error/.+ ",
                            paste0("-o ",hps_dir,"log_dir/log/",batch_name," -e ",hps_dir,"log_dir/error/",batch_name," -J ",batch_name," $SHELL "),
                            bsub_temp))
    hps_sh_file = paste0(LSF_project_path, sh_file)
    bsub_temp = paste0(bsub_temp, hps_sh_file)
  })

  commands$run = bsub
  commands
}
