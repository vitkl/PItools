##' find fold_enrichment distribution (in any domain or other feature) for each protein in the network
##' @name foldEnrichmentDist
##' @description \code{foldEnrichmentDist} finds fold_enrichment distribution (in any domain or other feature) for each protein in the first column of net. For more details: \code{\link[MItools]{foldEnrichment}}
##' @param net data.table that specifies the PPI network:
##' 1. first column specifies a protein ("IDs_interactor_viral")
##' 2. second column specifies it's interacting partner which annotations are of interest ("IDs_interactor_human")
##' 3. third column specifies the degree of a protein in the first column ("IDs_interactor_viral_degree")
##' @param protein_annot data.table that specifies annotations for a protein in the second column of \code{net}
##' 1. first column specifies a protein ("IDs_interactor_human")
##' 2. second column specifies features of that protein ("IDs_domain_human")
##' 3. third column specifies background frequency of those features ("domain_frequency)
##' @param N number of times to run permutation of PPI network
##' @param cores specify how many cores to use for parallel processing, default (NULL) is to detect all cores on the machine and use all minus one. When using LSF cluster you must specify the number of cores to use because \code{\link[BiocGenerics]{detectCores}} doen't know how much cores you have requested from LSF (with bsub -n) and detects all cores on the actual physical node.
##' @param seed seed for RNG for reproducible sampling
##' @param frequency fold enrichment or frequency in a set (if TRUE - frequency), also passed to \code{\link[MItools]{foldEnrichment}}
##' @return data.table containing fold enrichment for each domain (or other feature) - protein pair (2 columns: IDs_interactor_viral and sampled_fold_enrichment or sampled_domain_frequency_per_set)
##' @author Vitalii Kleshchevnikov
##' @import data.table
##' @import BiocGenerics
##' @usage
##' {
##' library(data.table)
##' data = fread("https://raw.githubusercontent.com/vitkl/viral_project/master/processed_data_files/viral_human_net_w_domains")
##' # generate minimal information tables
##' net = unique(data[,.(IDs_interactor_viral, IDs_interactor_human, IDs_interactor_viral_degree)])
##' protein_annot = unique(data[,.(IDs_interactor_human, IDs_domain_human, domain_frequency)])
##' viral_human_net_w_domains = unique(data[,.(IDs_interactor_viral, IDs_interactor_human, IDs_domain_human, domain_frequency_per_IDs_interactor_viral, fold_enrichment)])
##' foldEnrichmentDist(net, protein_annot, N = 10, cores = NULL, seed = 1, frequency = T)
##' }
foldEnrichmentDist = function(net, protein_annot, N = 1000, cores = NULL, seed = 1, frequency = T){

  # set up parallel processing
  # create cluster
  if(is.null(cores)) cores = detectCores()-1
  cl <- makeCluster(cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(data.table); library(MItools)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("net", "protein_annot", "frequency"), envir=environment())
  # set seed
  clusterSetRNGStream(cl, iseed = seed)
  if(!frequency){
    # run permutations N times
    fold_enrichment_dist = parReplicate(cl, n = N,
                                        expr = {
                                          net[, IDs_interactor_human := sample(IDs_interactor_human)]
                                          sample_net = foldEnrichment(net, protein_annot, frequency)

                                          # separate interactions involving no domains: otherwise these will get collapsed into one per viral protein
                                          sample_net_NA = sample_net[is.na(IDs_domain_human), .(IDs_interactor_viral, sample_fold_enrichment = fold_enrichment)]
                                          # first, collapse (unique()) counts by IDs_interactor_viral AND IDs_domain_human, then remove IDs_domain_human: otherwise identical values of counts will get collapsed into one (for example, when two domains have the same frequency among interactors of one viral protein)
                                          sample_net_notNA = unique(sample_net[!is.na(IDs_domain_human), .(IDs_interactor_viral, IDs_domain_human, sample_fold_enrichment = fold_enrichment)])[, IDs_domain_human := NULL]

                                          # rbind no domain and domain present cases, collapse all counts per viral protein into one cell
                                          fold_enrichment_dist = rbind(sample_net_NA, sample_net_notNA)[, .(sample_fold_enrichment = list(sample_fold_enrichment)), by = IDs_interactor_viral]

                                          setorder(fold_enrichment_dist, IDs_interactor_viral)
                                          fold_enrichment_dist_v = fold_enrichment_dist$sample_fold_enrichment
                                          names(fold_enrichment_dist_v) = fold_enrichment_dist$IDs_interactor_viral
                                          fold_enrichment_dist_v
                                        },
                                        simplify=TRUE, USE.NAMES=TRUE)
    fold_enrichment_dist = as.data.table(fold_enrichment_dist, keep.rownames = "IDs_interactor_viral")
    fold_enrichment_dist = fold_enrichment_dist[, .(sampled_fold_enrichment = as.numeric((unlist(.SD)))), .SDcols = paste0("V", 1:N), by = IDs_interactor_viral]
  }
  if(frequency){
    # run permutations N times
    fold_enrichment_dist = parReplicate(cl, n = N,
                                        expr = {
                                          net[, IDs_interactor_human := sample(IDs_interactor_human)]
                                          sample_net = foldEnrichment(net, protein_annot, frequency)

                                          # separate interactions involving no domains: otherwise these will get collapsed into one per viral protein
                                          sample_net_NA = sample_net[is.na(IDs_domain_human), .(IDs_interactor_viral, sample_domain_frequency_per_set = domain_frequency_per_set)]
                                          # first, collapse (unique()) counts by IDs_interactor_viral AND IDs_domain_human, then remove IDs_domain_human: otherwise identical values of counts will get collapsed into one (for example, when two domains have the same frequency among interactors of one viral protein)
                                          sample_net_notNA = unique(sample_net[!is.na(IDs_domain_human), .(IDs_interactor_viral, IDs_domain_human, sample_domain_frequency_per_set = domain_frequency_per_set)])[, IDs_domain_human := NULL]

                                          # rbind no domain and domain present cases, collapse all counts per viral protein into one cell
                                          fold_enrichment_dist = rbind(sample_net_NA, sample_net_notNA)[, .(sample_domain_frequency_per_set = list(sample_domain_frequency_per_set)), by = IDs_interactor_viral]

                                          setorder(fold_enrichment_dist, IDs_interactor_viral)
                                          fold_enrichment_dist_v = fold_enrichment_dist$sample_domain_frequency_per_set
                                          names(fold_enrichment_dist_v) = fold_enrichment_dist$IDs_interactor_viral
                                          fold_enrichment_dist_v
                                        },
                                        simplify=TRUE, USE.NAMES=TRUE)

    fold_enrichment_dist = as.data.table(fold_enrichment_dist, keep.rownames = "IDs_interactor_viral")
    fold_enrichment_dist = fold_enrichment_dist[, .(sampled_domain_frequency_per_set = as.numeric((unlist(.SD)))), .SDcols = paste0("V", 1:N), by = IDs_interactor_viral]
  }

  # stop the cluster
  stopCluster(cl)

  return(fold_enrichment_dist)
}

###########################################################################################

foldEnrichmentDist_pasteSplit = function(net, protein_annot, N = 1000, cores = NULL, seed = 1, frequency = T){

  # set up parallel processing
  # create cluster
  if(is.null(cores)) cores = detectCores()-1
  cl <- makeCluster(cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(data.table); library(MItools)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("net", "protein_annot", "frequency"), envir=environment())
  # set seed
  clusterSetRNGStream(cl, iseed = seed)
  if(!frequency){
    # run permutations N times
    fold_enrichment_dist = parReplicate(cl, n = N,
                                        expr = {
                                          net[, IDs_interactor_human := sample(IDs_interactor_human)]
                                          sample_net = foldEnrichment(net, protein_annot, frequency)

                                          # separate interactions involving no domains: otherwise these will get collapsed into one per viral protein
                                          sample_net_NA = sample_net[is.na(IDs_domain_human), .(IDs_interactor_viral, sample_fold_enrichment = fold_enrichment)]
                                          # first, collapse (unique()) counts by IDs_interactor_viral AND IDs_domain_human, then remove IDs_domain_human: otherwise identical values of counts will get collapsed into one (for example, when two domains have the same frequency among interactors of one viral protein)
                                          sample_net_notNA = unique(sample_net[!is.na(IDs_domain_human), .(IDs_interactor_viral, IDs_domain_human, sample_fold_enrichment = fold_enrichment)])[, IDs_domain_human := NULL]

                                          # rbind no domain and domain present cases, collapse all counts per viral protein into one cell
                                          fold_enrichment_dist = rbind(sample_net_NA, sample_net_notNA)[, .(sample_fold_enrichment = paste0(sample_fold_enrichment, collapse = "|")), by = IDs_interactor_viral]

                                          setorder(fold_enrichment_dist, IDs_interactor_viral)
                                          fold_enrichment_dist_v = fold_enrichment_dist$sample_fold_enrichment
                                          names(fold_enrichment_dist_v) = fold_enrichment_dist$IDs_interactor_viral
                                          fold_enrichment_dist_v
                                        },
                                        simplify=TRUE, USE.NAMES=TRUE)
    fold_enrichment_dist = as.data.table(fold_enrichment_dist, keep.rownames = "IDs_interactor_viral")
    fold_enrichment_dist = fold_enrichment_dist[, .(sampled_fold_enrichment = as.numeric(unlist(strsplit(unlist(.SD), "\\|")))), .SDcols = paste0("V", 1:N), by = IDs_interactor_viral]
  }
  if(frequency){
    # run permutations N times
    fold_enrichment_dist = parReplicate(cl, n = N,
                                        expr = {
                                          net[, IDs_interactor_human := sample(IDs_interactor_human)]
                                          sample_net = foldEnrichment(net, protein_annot, frequency)

                                          # separate interactions involving no domains: otherwise these will get collapsed into one per viral protein
                                          sample_net_NA = sample_net[is.na(IDs_domain_human), .(IDs_interactor_viral, sample_domain_frequency_per_set = domain_frequency_per_set)]
                                          # first, collapse (unique()) counts by IDs_interactor_viral AND IDs_domain_human, then remove IDs_domain_human: otherwise identical values of counts will get collapsed into one (for example, when two domains have the same frequency among interactors of one viral protein)
                                          sample_net_notNA = unique(sample_net[!is.na(IDs_domain_human), .(IDs_interactor_viral, IDs_domain_human, sample_domain_frequency_per_set = domain_frequency_per_set)])[, IDs_domain_human := NULL]

                                          # rbind no domain and domain present cases, collapse all counts per viral protein into one cell
                                          fold_enrichment_dist = rbind(sample_net_NA, sample_net_notNA)[, .(sample_domain_frequency_per_set = paste0(sample_domain_frequency_per_set, collapse = "|")), by = IDs_interactor_viral]

                                          setorder(fold_enrichment_dist, IDs_interactor_viral)
                                          fold_enrichment_dist_v = fold_enrichment_dist$sample_domain_frequency_per_set
                                          names(fold_enrichment_dist_v) = fold_enrichment_dist$IDs_interactor_viral
                                          fold_enrichment_dist_v
                                        },
                                        simplify=TRUE, USE.NAMES=TRUE)

    fold_enrichment_dist = as.data.table(fold_enrichment_dist, keep.rownames = "IDs_interactor_viral")
    fold_enrichment_dist = fold_enrichment_dist[, .(sampled_domain_frequency_per_set = as.numeric(unlist(strsplit(unlist(.SD), "\\|")))), .SDcols = paste0("V", 1:N), by = IDs_interactor_viral]
  }

  # stop the cluster
  stopCluster(cl)

  return(fold_enrichment_dist)
}
