##' Generate random GRanges of equal length in the same sequence
##' @rdname randomGRanges
##' @name randomGRanges
##' @author Vitalii Kleshchevnikov
##' @param GRanges GRanges object
##' @param N how many randomised GRanges to create
##' @param replace passed to \code{\link[base]{sample.int}}
##' @param within1sequence logical, if TRUE ranges of equal lenght are sampled from the same sequence. If false, ranges of equal length are sampled from all sequences (while keeping the number of ranges in each sequence).
##' @return Genomic Ranges list object containing \code{N} randomised GRanges objects
##' @import GenomicRanges
##' @export randomGRanges
##' @seealso \code{\link{ELMdb2GRanges}}, \code{\link{GRangesINinteractionSubsetFASTA}}
randomGRanges = function(GRanges, N = 1, replace = T, within1sequence = T){
  widths = width(GRanges)
  seq_names = as.character(seqnames(GRanges))
  seq_lengths = seqlengths(GRanges)
  GRanges_length = length(GRanges)

  if(within1sequence){
    max_starts = seq_lengths[seq_names] - widths
    starts = sapply(max_starts, function(max_start) {
      sample.int(n = max_start, size = N, replace = replace)
    })
    if(N == 1) starts = matrix(starts, nrow = 1, ncol = GRanges_length)

    res = apply(starts, 1, function(start){
      temp_res = GRanges(seqnames=seq_names,
              ranges=IRanges(start=start, end=NULL, width=widths, names=seq_names),
              strand=NULL,
              seqlengths = seq_lengths)
      temp_res$source = paste0("random_within1sequence_", within1sequence)
      temp_res$type = "sequence_feature"
      temp_res$for_benchmarking = 0
      colnames_ = colnames(mcols(GRanges))
      colnames_ = colnames_[!colnames_ %in% c("source", "type", "for_benchmarking")]
      for(colname in colnames_){
        eval(parse(text = paste0("temp_res$",colname," = NA")))
      }
      mcols(temp_res) = mcols(temp_res)[,colnames(mcols(GRanges))]
      temp_res
    })
  } else {
    res = replicate(N, {
      seq_names = seq_names[sample.int(length(seq_names))]
      max_starts = seq_lengths[seq_names] - widths
      start = sapply(max_starts, function(max_start) {
        sample.int(n = max_start, size = 1, replace = replace)
      })
      temp_res = GRanges(seqnames=seq_names,
              ranges=IRanges(start=start, end=NULL, width=widths, names=seq_names),
              strand=NULL,
              seqlengths = seq_lengths)
      temp_res$source = paste0("random_within1sequence_", within1sequence)
      temp_res$type = "sequence_feature"
      temp_res$for_benchmarking = 0
      colnames_ = colnames(mcols(GRanges))
      colnames_ = colnames_[!colnames_ %in% c("source", "type", "for_benchmarking")]
      for(colname in colnames_){
        eval(parse(text = paste0("temp_res$",colname," = NA")))
      }
      mcols(temp_res) = mcols(temp_res)[,colnames(mcols(GRanges))]
      temp_res
    })
  }
  return(res)
}
