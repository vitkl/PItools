##' Find occuring GRanges overlapping with the benchmarking set (positive + negative examples)
##' @rdname findOverlapsBench
##' @name findOverlapsBench
##' @author Vitalii Kleshchevnikov
##' @param occuring GRanges object, occuring feature instances
##' @param benchmarking GRanges object, feature instances for benchmarking
##' @param predictor_col name of the metadata column in \code{occuring} that contains 1 - predictor value (such as p-value)
##' @param labels_col name of the metadata column in \code{benchmarking} that contains labels for benchmarking (for, example 0 for negative and 1 for positive)
##' @param normalise logical, normalise predictor value, just in case predictor doesn't span the full range between 0 ... 1
##' @param ... arguments such as \code{maxgap}, \code{minoverlap} passed to \code{\link[GenomicRanges]{findOverlaps}}
##' @return Genomic Ranges list object containing \code{N} randomised GRanges objects
##' @import GenomicRanges
##' @export findOverlapsBench
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{GRangesINinteractionSubsetFASTA}}
findOverlapsBench = function(occuring, benchmarking, predictor_col = "Sig", labels_col = "for_benchmarking", normalise = T, select_predictor_per_range = min, ...){
  suppressWarnings({ #suppress warnings about different sequences being present in the occuring and benchmarking sets
    overlap_inds = findOverlaps(occuring, benchmarking, ...)
  })
  occuring_subset = occuring[queryHits(overlap_inds)]
  # benchmarkings_subset = benchmarkings[subjectHits(overlap_inds)]
  mcols(occuring_subset) = c(mcols(occuring_subset), mcols(benchmarking)[subjectHits(overlap_inds),])

  suppressWarnings({
    no_overlap = subsetByOverlaps(benchmarking, occuring, invert = T, ...)
  })

  # normalise, just in case predictor doesn't span the full range 0 ... 1
  if(normalise) mcols(occuring_subset)[,predictor_col] = mcols(occuring_subset)[,predictor_col]/max(mcols(occuring_subset)[,predictor_col])

  res_overlap = data.frame(predictions = 1 - mcols(occuring_subset)[,predictor_col],
                           labels = mcols(occuring_subset)[,labels_col],
                           stringsAsFactors = F)
  if(length(no_overlap) >= 1){
    res_no_overlap = data.frame(predictions = 0,
                                labels = mcols(no_overlap)[,labels_col],
                                stringsAsFactors = F)
  } else {
    res_no_overlap = data.frame(predictions = numeric(),
                                labels = numeric(),
                                stringsAsFactors = F)
  }

  # choose only one p-value for the same range: occuring_subset
  res_overlap_list = split(res_overlap, paste0(seqnames(occuring_subset), start(occuring_subset), end(occuring_subset)))
  res_overlap_list = lapply(res_overlap_list, function(res_over) {
    res_over$predictions = select_predictor_per_range(res_over$predictions)
    unique(res_over)
  })
  res_overlap = Reduce(rbind, res_overlap_list)
  # choose only one p-value for the same range: no_overlap - useless code: predictions are 0
  #res_no_overlap_list = split(res_no_overlap, paste0(seqnames(no_overlap), start(no_overlap), end(no_overlap)))
  #res_no_overlap_list = lapply(res_no_overlap_list, function(res_over) {
  #  res_over$predictions = select_predictor_per_range(res_over$predictions)
  #  unique(res_over)
  #})
  #res_no_overlap = Reduce(rbind, res_no_overlap_list)

  res = rbind(res_overlap, res_no_overlap)
  list(overlapping_GRanges = occuring_subset, for_ROC = res, not_found = no_overlap)
}
