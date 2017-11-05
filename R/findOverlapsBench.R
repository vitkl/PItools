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
findOverlapsBench = function(occuring, benchmarking, predictor_col = "Sig", labels_col = "for_benchmarking", normalise = T, ...){
  overlap_inds = findOverlaps(occuring, benchmarking, ...)
  occuring_subset = occuring[queryHits(overlap_inds)]
  # benchmarkings_subset = benchmarkings[subjectHits(overlap_inds)]
  mcols(occuring_subset) = c(mcols(occuring_subset), mcols(benchmarking)[subjectHits(overlap_inds),])

  no_overlap = subsetByOverlaps(benchmarking, occuring, invert = T, ...)

  # normalise, just in case predictor doesn't span the full range 0 ... 1
  if(normalise) mcols(occuring_subset)[,predictor_col] = mcols(occuring_subset)[,predictor_col]/max(mcols(occuring_subset)[,predictor_col])

  res_overlap = data.frame(predictions = 1 - mcols(occuring_subset)[,predictor_col], labels = mcols(occuring_subset)[,labels_col])
  res_no_overlap = data.frame(predictions = 0, labels = mcols(no_overlap)[,labels_col])
  res = rbind(res_overlap, res_no_overlap)
  list(overlapping_GRanges = occuring_subset, for_ROC = res)
}
