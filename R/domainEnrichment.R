##' domainEnrichment: filter domains by background domain count => foldEnrichmentDist => foldEnrichmentPval => pval correction => plot pvals
##' @param backgr_domain_count discard domains that have background count of less or equal \code{backgr_domain_count}
##' @param all.data data.table: all data (net + protein_annot + data) + additional columns, domain_count is necessary for the function to work
##' @param net data.table that specifies the PPI network:
##' 1. first column specifies a protein ("IDs_interactor_viral")
##' 2. second column specifies it's interacting partner which annotations are of interest ("IDs_interactor_human")
##' 3. third column specifies the degree of a protein in the first column ("IDs_interactor_viral_degree")
##' @param protein_annot data.table that specifies annotations for a protein in the second column of \code{net}
##' 1. first column specifies a protein ("IDs_interactor_human")
##' 2. second column specifies features of that protein ("IDs_domain_human")
##' 3. third column specifies background frequency of those features ("domain_frequency)
##' @param frequency fold enrichment or frequency in a set (if TRUE - frequency)
##' @param N number of times to run permutation of PPI network
##' @param cores specify how many cores to use for parallel processing, default (NULL) is to detect all cores on the machine and use all minus one. When using LSF cluster you must specify the number of cores to use because \code{\link[BiocGenerics]{detectCores}} doen't know how much cores you have requested from LSF (with bsub -n) and detects all cores on the actual physical node.
##' @param seed seed for RNG for reproducible sampling
##' @import data.table
##' @import qvalue
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_histogram
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 theme_light
##' @importFrom ggplot2 xlim
##' @export domainEnrichment
domainEnrichment = function(backgr_domain_count, all.data, net, protein_annot, data, N = 1000, cores = NULL, seed = 1, frequency = T){
  # filter out domains with lower background domain count than backgr_domain_count
  domains2keep = all.data[domain_count > backgr_domain_count, IDs_domain_human]

  protein_annot = protein_annot[IDs_domain_human %in% domains2keep,]
  data = data[IDs_domain_human %in% domains2keep,]
  all.data = all.data[IDs_domain_human %in% domains2keep,]

  viral_foldEnrichDist = foldEnrichmentDist(net = net,
                                            protein_annot = protein_annot,
                                            N = N, cores = cores, seed = seed)

  Pvals = foldEnrichmentPval(fold_enrichment_dist = viral_foldEnrichDist,
                             data = data, cores = cores)

  Pvals[, Pval_fdr := p.adjust(Pval, method = "fdr")]
  Pvals[, Qval := qvalue(Pval)$qvalues]

  ggplot(Pvals, aes(x = Pval_fdr)) + geom_histogram(bins = 100) + ggtitle(paste0("viral protein and human domain association \n FDR adjusted pvalue distribution \n background domain count > ", backgr_domain_count)) + theme_light() + xlim(0,1)

  # merge results to original data
  # all.data = all.data[Pvals, on = c("IDs_interactor_viral", "IDs_domain_human", "fold_enrichment")]

  return(list(Pvals,all.data))
}
