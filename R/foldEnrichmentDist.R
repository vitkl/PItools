##' \code{foldEnrichmentDist} finds fold_enrichment distribution (in any domain or other feature) for each protein
##' @name foldEnrichmentDist
##' @description \code{foldEnrichmentDist} finds fold_enrichment distribution (in any domain or other feature) for each protein in the first column of net. For more details: \code{\link{foldEnrichment}}
##' @param net data.table that specifies the PPI network:
##' 1. first column specifies a protein ("IDs_interactor_viral")
##' 2. second column specifies it's interacting partner which annotations are of interest ("IDs_interactor_human")
##' 3. third column specifies the degree of a protein in the first column ("IDs_interactor_viral_degree")
##' @param protein_annot data.table that specifies annotations for a protein in the second column of \code{net}
##' 1. first column specifies a protein ("IDs_interactor_human")
##' 2. second column specifies features of that protein ("IDs_domain_human")
##' 3. third column specifies background frequency of those features ("domain_frequency)
##' @return data.table containing fold enrichment for each domain - protein pair
##' @author Vitalii Kleshchevnikov
##' @import data.table
##' @import BiocGenerics
##' @export foldEnrichmentDist
foldEnrichmentDist = function(net, protein_annot, N = 1000, cores = NULL){

  # set up parallel processing
  # create cluster
  if(is.null(cores)) cores = detectCores()-1
  cl <- makeCluster(cores)
  # get library support needed to run the code
  clusterEvalQ(cl, {library(data.table); library(MItools)})
  # put objects in place that might be needed for the code
  clusterExport(cl, c("net", "protein_annot"), envir=environment())
  # setorder(net, IDs_interactor_viral)
  fold_enrichment_dist = parReplicate(cl, n = N,
               expr = {
                 net[, IDs_interactor_human := sample(IDs_interactor_human)]
                 sample_net = foldEnrichment(net, protein_annot)
                 # find fold_enrichment distribution in any domain for each viral protein
                 fold_enrichment_dist = unique(sample_net[,.(IDs_interactor_viral, sample_fold_enrichment = fold_enrichment)])[, .(sample_fold_enrichment = paste0(sample_fold_enrichment, collapse = "|")), by = IDs_interactor_viral]
                 setorder(fold_enrichment_dist, IDs_interactor_viral)
                 fold_enrichment_dist_v = fold_enrichment_dist$sample_fold_enrichment
                 names(fold_enrichment_dist_v) = fold_enrichment_dist$IDs_interactor_viral
                 fold_enrichment_dist_v
               },
               simplify=TRUE, USE.NAMES=TRUE)

  # stop the cluster
  stopCluster(cl)

  fold_enrichment_dist = as.data.table(fold_enrichment_dist, keep.rownames = "IDs_interactor_viral")
  fold_enrichment_dist = fold_enrichment_dist[, .(sampled_fold_enrichment = as.numeric(unlist(strsplit(unlist(.SD), "\\|")))), .SDcols = paste0("V", 1:N), by = IDs_interactor_viral]

  return(fold_enrichment_dist)
}
