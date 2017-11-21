##' Find linear motifs (QSLIMFinder or SLIMFinder) in the protein interaction network
##' @rdname PPInetwork2SLIMFinder
##' @name PPInetwork2SLIMFinder
##' @author Vitalii Kleshchevnikov
##' @param dataset_name refer to \code{\link[MItools]{mBenchmarkMotifs}}
##' @param interaction_main_set clean_MItab class, use this set of protein interactions to construct QSLIMFinder datasets
##' @param interaction_query_set  clean_MItab class, use this set of protein interactions as a query (+ add to the QSLIMFinder datasets). SLIMFinder \code{analysis_type} also requires this option because it add proteins from these interactions to the SLIMFinder datasets
##' @param analysis_type "qslimfinder" or "slimfinder"
##' @param options any options from QSLIMFinder or SLIMFinder. Detail http://rest.slimsuite.unsw.edu.au/docs&page=module:qslimfinder or http://rest.slimsuite.unsw.edu.au/docs&page=module:slimfinder => Commandline
##' @param path2domain_enrich relative path to domain enrichment results RData
##' @param domain_enrich_object which object contains domain enrichment results in \code{path2domain_enrich}?
##' @param seed4datasets_col column in domain_enrich_object that specified seed for QSLIMFinder datasets
##' @param fasta_path relative path (from the project folder) to the FASTA file containing sequences for all proteins in \code{interaction_main_set} and \code{interaction_query_set}
##' @param main_set_only logical, sequence set for motif search should contain only proteins from \code{interaction_main_set}. If FALSE, non-query proteins from \code{interaction_query_set} are also included. Argument for \code{\link{listInteractionSubsetFASTA}}
##' @param domain_pvalue_cutoff construct SLIMFinder datasets using interactions of proteins that contain domain associated to protein in the query set with p-value \code{domain_pvalue_cutoff} or lower
##' @param SLIMFinder_dir directory to store SLIMFinder datasets and results within the project directory
##' @param LSF_project_path full path to the project directory
##' @param software_path relative path (from the project folder) to the directory containing slimsuite, blast, iupred # "../software/cluster/" or "../software/"
##' @param length_set1_min mininal number of proteins in a QSLIMFinder dataset from \code{interaction_main_set}. Argument for \code{\link{filterInteractionSubsetFASTA_list}}
##' @param length_set2_min mininal number of proteins in a QSLIMFinder dataset from \code{interaction_query_set}. Argument for \code{\link{filterInteractionSubsetFASTA_list}}
##' @param write_log FALSE will not allow runQSLIMFinder to detect crashed jobs
##' @param N_seq number of sequences per batch
##' @details QSLIMFinder command line options (http://rest.slimsuite.unsw.edu.au/docs&page=module:qslimfinder)
##'
##'### Basic Input/Output Options ###
##'
##'seqin=FILE : Sequence file to search [None]
##'
##'batch=LIST : List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
##'
##'query=LIST : Return only SLiMs that occur in 1+ Query sequences (Name/AccNum/Seq Number) [1]
##'
##'addquery=FILE : Adds query sequence(s) to batch jobs from FILE [None]
##'
##'maxseq=X : Maximum number of sequences to process [500]
##'
##'maxupc=X : Maximum UPC size of dataset to process [0]
##'
##'sizesort=X : Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
##'
##'walltime=X : Time in hours before program will abort search and exit [1.0]
##'
##'resfile=FILE : Main QSLiMFinder results table [qslimfinder.csv]
##'
##'resdir=PATH : Redirect individual output files to specified directory (and look for intermediates) [QSLiMFinder/]
##'
##'buildpath=PATH : Alternative path to look for existing intermediate files [SLiMFinder/]
##'
##'force=T/F : Force re-running of BLAST, UPC generation and SLiMBuild [False]
##'
##'pickup=T/F : Pick-up from aborted batch run by identifying datasets in resfile using RunID [False]
##'
##'dna=T/F : Whether the sequences files are DNA rather than protein [False]
##'
##'alphabet=LIST : List of characters to include in search (e.g. AAs or NTs) [default AA or NT codes]
##'
##'megaslim=FILE : Make/use precomputed results for a proteome (FILE) in fasta format [None]
##'
##'megablam=T/F : Whether to create and use all-by-all GABLAM results for (gablamdis) UPC generation [False]
##'
##'ptmlist=LIST : List of PTM letters to add to alphabet for analysis and restrict PTM data []
##'
##'ptmdata=DSVFILE : File containing PTM data, including AccNum, ModType, ModPos, ModAA, ModCode
##'
##'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##'
##'SLiMBuild Options I
##'
##'efilter=T/F : Whether to use evolutionary filter [True]
##'
##'blastf=T/F : Use BLAST Complexity filter when determining relationships [True]
##'
##'blaste=X : BLAST e-value threshold for determining relationships [1e=4]
##'
##'altdis=FILE : Alternative all by all distance matrix for relationships [None]
##'
##'gablamdis=FILE : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
##'
##'homcut=X : Max number of homologues to allow (to reduce large multi-domain families) [0]
##'
##'SLiMBuild Options II
##'
##'masking=T/F : Master control switch to turn off all masking if False [True]
##'
##'dismask=T/F : Whether to mask ordered regions (see rje_disorder for options) [False]
##'
##'consmask=T/F : Whether to use relative conservation masking [False]
##'
##'ftmask=LIST : UniProt features to mask out (True=EM,DOMAIN,TRANSMEM) []
##'
##'imask=LIST : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
##'
##'compmask=X,Y : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
##'
##'casemask=X : Mask Upper or Lower case [None]
##'
##'motifmask=X : List (or file) of motifs to mask from input sequences []
##'
##'metmask=T/F : Masks the N-terminal M (can be useful if termini=T) [True]
##'
##'posmask=LIST : Masks list of position-specific aas, where list = pos1:aas,pos2:aas [2:A]
##'
##'aamask=LIST : Masks list of AAs from all sequences (reduces alphabet) []
##'
##'qregion=X,Y : Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
##'
##'SLiMBuild Options III
##'
##'termini=T/F : Whether to add termini characters (^ & $) to search sequences [True]
##'
##'minwild=X : Minimum number of consecutive wildcard positions to allow [0]
##'
##'maxwild=X : Maximum number of consecutive wildcard positions to allow [2]
##'
##'slimlen=X : Maximum length of SLiMs to return (no. non-wildcard positions) [5]
##'
##'minocc=X : Minimum number of unrelated occurrences for returned SLiMs. (Proportion of UP if < 1) [0.05]
##'
##'absmin=X : Used if minocc<1 to define absolute min. UP occ [3]
##'
##'alphahelix=T/F : Special i, i+3/4, i+7 motif discovery [False]
##'
##'SLiMBuild Options IV
##'
##'ambiguity=T/F : (preamb=T/F) Whether to search for ambiguous motifs during motif discovery [True]
##'
##'ambocc=X : Min. UP occurrence for subvariants of ambiguous motifs (minocc if 0 or > minocc) [0.05]
##'
##'absminamb=X : Used if ambocc<1 to define absolute min. UP occ [2]
##'
##'equiv=LIST : List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,FYH,KRH,DE,ST]
##'
##'wildvar=T/F : Whether to allow variable length wildcards [True]
##'
##'combamb=T/F : Whether to search for combined amino acid degeneracy and variable wildcards [False]
##'
##'SLiMBuild Options V
##'
##'musthave=LIST : Returned motifs must contain one or more of the AAs in LIST (reduces search space) []
##'
##'focus=FILE : FILE containing focal groups for SLiM return (see Manual for details) [None]
##'
##'focusocc=X : Motif must appear in X+ focus groups (0 = all) [0]
##'
##'* See also rje_slimcalc options for occurrence-based calculations and filtering *
##'
##'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##'
##'### SLiMChance Options ###
##'
##'cloudfix=T/F : Restrict output to clouds with 1+ fixed motif (recommended) [False]
##'
##'slimchance=T/F : Execute main QSLiMFinder probability method and outputs [True]
##'
##'sigprime=T/F : Calculate more precise (but more computationally intensive) statistical model [False]
##'
##'sigv=T/F : Use the more precise (but more computationally intensive) fix to mean UPC probability [False]
##'
##'qexact=T/F : Calculate exact Query motif space (True) or over-estimate from dimers (False) (quicker) [True]
##'
##'probcut=X : Probability cut-off for returned motifs [0.1]
##'
##'maskfreq=T/F : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [False]
##'
##'aafreq=FILE : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
##'
##'aadimerfreq=FILE: Use empirical dimer frequencies from FILE (fasta or *.aadimer.tdt) (!!!Experimental!!!) [None]
##'
##'negatives=FILE : Multiply raw probabilities by under-representation in FILE (!!!Experimental!!!) [None]
##'
##'smearfreq=T/F : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
##'
##'seqocc=T/F : Whether to upweight for multiple occurrences in same sequence (heuristic) [False]
##'
##'probscore=X : Score to be used for probability cut-off and ranking (Prob/Sig) [Sig]
##'
##'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##'
##'Advanced Output Options I
##'
##'clouds=X : Identifies motif "clouds" which overlap at 2+ positions in X+ sequences (0=minocc / -1=off) [2]
##'
##'runid=X : Run ID for resfile (allows multiple runs on same data) [DATE:TIME]
##'
##'logmask=T/F : Whether to log the masking of individual sequences [True]
##'
##'slimcheck=FILE : Motif file/list to add to resfile output []
##'
##'Advanced Output Options II
##'
##'teiresias=T/F : Replace TEIRESIAS, making *.out and *.mask.fasta files [False]
##'
##'slimdisc=T/F : Emulate SLiMDisc output format (*.rank & *.dat.rank + TEIRESIAS *.out & *.fasta) [False]
##'
##'extras=X : Whether to generate additional output files (alignments etc.) [1]
##'
##'--1 = No output beyond main results file
##'
##'- 0 = Generate occurrence file and cloud file
##'
##'- 1 = Generate occurrence file, alignments and cloud file
##'
##'- 2 = Generate all additional QSLiMFinder outputs
##'
##'- 3 = Generate SLiMDisc emulation too (equiv extras=2 slimdisc=T)
##'
##'targz=T/F : Whether to tar and zip dataset result files (UNIX only) [False]
##'
##'savespace=0 : Delete "unneccessary" files following run (best used with targz): [0]
##'
##'- 0 = Delete no files
##'
##'- 1 = Delete all bar *.upc and *.pickle
##'
##'- 2 = Delete all bar *.upc (pickle added to tar)
##'
##'- 3 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)
##'
##'Advanced Output Options III
##'
##'topranks=X : Will only output top X motifs meeting probcut [1000]
##'
##'minic=X : Minimum information content for returned motifs [2.1]
##'
##'allsig=T/F : Whether to also output all SLiMChance combinations (Sig/SigV/SigPrime/SigPrimeV) [False]
##'
##'
##' @return path to RData containing all objects used by this pipeline
##' @import data.table
##' @import rtracklayer
##' @import Biostrings
##' @export PPInetwork2SLIMFinder
##' @seealso \code{\link{listInteractionSubsetFASTA}}, \code{\link{runQSLIMFinder}}, \code{\link{groupQSLIMFinderCommand}}, \code{\link{mQSLIMFinderCommand}}, \code{\link{runCompariMotif3}}, \code{\link{readQSLIMFinderMain}}, \code{\link{readQSLIMFinderOccurence}}, \code{\link{writeInteractionSubsetFASTA_list}}, \code{\link{domainProteinPairMatch}}, \code{\link{filterInteractionSubsetFASTA_list}}, \code{\link{removeInteractionNoFASTA}}
PPInetwork2SLIMFinder = function(dataset_name = "SLIMFinder",
                                 interaction_main_set = all_human_interaction,
                                 interaction_query_set = all_viral_interaction,
                                 analysis_type = "qslimfinder",
                                 options = "dismask=T consmask=F cloudfix=T probcut=0.3 minwild=0 maxwild=2 slimlen=5 alphahelix=F maxseq=1500 savespace=1 iuchdir=T",
                                 path2domain_enrich = "./processed_data_files/what_we_find_VS_ELM_clust20171019.RData",
                                 domain_enrich_object = "res_count",
                                 seed4datasets_col = "IDs_interactor_human",
                                 fasta_path = "./data_files/all_human_viral_proteins.fasta",
                                 main_set_only = F,
                                 domain_pvalue_cutoff = 1,
                                 SLIMFinder_dir = paste0("./",dataset_name,"/"),
                                 LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                                 software_path = "../software/cluster/",
                                 length_set1_min = 2,
                                 length_set2_min = 1,
                                 write_log = T,
                                 N_seq = 200)
{
  # check class correctness
  if(!grepl("clean_MItab",class(MITABdata))) stop("MITABdata is not of class clean_MItab27 or related clean_MItab class")

  # load the domain analysis results
  load(path2domain_enrich, envir = environment())
  eval(parse(text = paste0("domain_res = ",domain_enrich_object)))
  rm(list = ls()[!ls() %in% c("domain_res", "dataset_name", "interaction_main_set",  "interaction_query_set", "analysis_type", "options", "path2domain_enrich", "domain_enrich_object",  "seed4datasets_col", "fasta_path", "main_set_only", "domain_pvalue_cutoff", "SLIMFinder_dir", "LSF_project_path", "software_path", "length_set1_min", "length_set2_min", "write_log", "N_seq")], envir = environment())

  # choose pvalue cutoff:
  plot(domain_res)
  eval(parse(text = paste0("proteins_w_signif_domains = unique(domain_res$data_with_pval[p.value <= domain_pvalue_cutoff, ", seed4datasets_col,"])")))

  # load FASTA
  all.fasta = readAAStringSet(filepath = fasta_path, format = "fasta")

  # remove interactions if one of the interactors does not have a FASTA sequence
  interaction_main_set = removeInteractionNoFASTA(interaction_main_set, all.fasta)
  interaction_query_set = removeInteractionNoFASTA(interaction_query_set, all.fasta)

  forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = interaction_main_set,
                                             interaction_set2 = interaction_query_set,
                                             seed_id_vect = proteins_w_signif_domains,
                                             fasta = all.fasta,
                                             single_interact_from_set2 = T, set1_only = main_set_only)

  forSLIMFinder_Ready = filterInteractionSubsetFASTA_list(forSLIMFinder,  length_set1_min = length_set1_min, length_set2_min = length_set2_min)

  # filter for only sets where seed protein - query protein pair matches significant domain - query protein pair
  domain_filt = domain_res
  domain_filt$data_with_pval = domain_filt$data_with_pval[p.value <= domain_pvalue_cutoff,]
  forSLIMFinder_Ready = domainProteinPairMatch(forSLIMFinder_Ready, domain_filt, remove = T)

  # write datasets (fasta + query)
  if(!dir.exists(SLIMFinder_dir)) dir.create(SLIMFinder_dir)
  forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
                                                             dir = SLIMFinder_dir, analysis_type = analysis_type)

  # create bash commands that will run QSLIMFinder
  all_commands = mQSLIMFinderCommand(file_list = forSLIMFinder_file_list,
                                     slimpath = paste0(software_path, "slimsuite/tools/"),
                                     blast = paste0(software_path, "ncbi_blast_2.6.0/bin/"),
                                     iupred = paste0(software_path, "iupred/iupred"),
                                     options = options,
                                     LSF_cluster_mode = T,
                                     LSF_project_path = LSF_project_path,
                                     LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                                     analysis_type = analysis_type,
                                     write_log = write_log)

  # group commands so that each run contains large number of sequences
  all_commands = groupQSLIMFinderCommand(commands = all_commands,
                                         InteractionSubsetFASTA_list = forSLIMFinder_Ready,
                                         sh_dir = "/sh_dir/",
                                         LSF_project_path = LSF_project_path,
                                         dataset_name = dataset_name, N_seq = N_seq, write_log = write_log)

  runQSLIMFinder(commands = all_commands, file_list = forSLIMFinder_file_list, onLSF = T)

  # read and bring together results
  resultdir = paste0(SLIMFinder_dir, "result/")
  if(!dir.exists(resultdir)) dir.create(resultdir)
  QSLIMFinder_main_result = readQSLIMFinderMain(outputfile = forSLIMFinder_file_list$outputfile)
  QSLIMFinder_occurence = readQSLIMFinderOccurence(outputdir = forSLIMFinder_file_list$outputdir)
  fwrite(QSLIMFinder_main_result, paste0(resultdir, "main_result.txt"), sep = "\t")

  fwrite(QSLIMFinder_occurence, paste0(resultdir, "occurence.txt"), sep = "\t")
  writePatternList(QSLIMFinder_main_result, filename = paste0(resultdir, "motifs.txt"))

  # compare motifs to ELM and to each other
  runCompariMotif3(input_file = paste0(resultdir, "motifs.txt"),
                   slimpath = paste0(software_path, "slimsuite/tools/"),
                   run = T, with = "db",
                   out_file = paste0(resultdir, "comparimotif.tdt"))
  runCompariMotif3(input_file = paste0(resultdir, "motifs.txt"),
                   slimpath = paste0(software_path, "slimsuite/tools/"),
                   run = T, with = "self",
                   out_file = paste0(resultdir, "comparimotif_with_self.tdt"))
  R.utils::gzip(paste0(resultdir, "comparimotif_with_self.compare.tdt"), paste0(resultdir, "comparimotif_with_self.compare.tdt.gz"))
  R.utils::gzip(paste0(resultdir, "comparimotif_with_self.compare.xgmml"), paste0(resultdir, "comparimotif_with_self.compare.xgmml.gz"))

  # compress input, remove input, log and output
  tar(paste0(SLIMFinder_dir, "input.gz"),paste0(SLIMFinder_dir, "input/"), compression='gzip')
  #tar(paste0(SLIMFinder_dir, "log_dir.gz"),paste0(SLIMFinder_dir, "log_dir/"), compression='gzip')
  #tar(paste0(SLIMFinder_dir, "output.gz"),paste0(SLIMFinder_dir, "output"), compression='gzip')
  unlink(paste0(SLIMFinder_dir, "input/"), recursive = T)
  unlink(paste0(SLIMFinder_dir, "log_dir/"))
  unlink(paste0(SLIMFinder_dir, "output/"))

  # save R session to RData
  AnalysisDate = Sys.Date()
  AnalysisSessionInfo = sessionInfo()

  filename = paste0("./processed_data_files/QSLIMFinder_instances_h2v_",dataset_name,"_clust",format(Sys.Date(), "%Y%m"),".RData")
  save(list = ls(), file=filename, envir = environment())
  return(filename)
}
