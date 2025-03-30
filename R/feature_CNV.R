#' @title Get CNV
#' @description
#' A function, is used to get CNV.
#'
#' @param path,bam bamfile which you want to read
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand
#' @importFrom GenomicAlignments cigar
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BSgenome getSeq
#' @importFrom tidyr separate
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract_all
#' @importFrom mgsub mgsub
#' @importFrom dplyr "%>%"
#' @return CNV
#' @export
#'
#'
feature_CNV <- function(path,bam){

  file <- paste0(path ,"/",bam)
  print(file)
  flag1 <- scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                       hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                       isFirstMateRead = NA, isSecondMateRead = NA,
                       isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                       isDuplicate = NA, isSupplementaryAlignment = NA)
  param1 <- ScanBamParam(
    flag=flag1,
    what=c("pos","qwidth"),
    tag = c("AS")
  )
  bam_data <- readGAlignments(file, use.names=TRUE, param=param1)


  regions <- list()
  current_region <- NULL
  all_reads_len = 0

  for (i in seq_along(bam_data)) {
    start_pos <- start(bam_data[i])
    end_pos <- end(bam_data[i])
    all_reads_len = all_reads_len + bam_data[i]@elementMetadata@listData$AS

    if (is.null(current_region)) {

      current_region <- c(start_pos, end_pos)
    } else if (start_pos <= current_region[2]) {

      current_region[2] <- max(current_region[2], end_pos)
    } else {

      regions <- append(regions, list(current_region))
      current_region <- c(start_pos, end_pos)
    }
  }


  if (!is.null(current_region)) {
    regions <- append(regions, list(current_region))
  }

  all_region_len = 0
  for (region in regions) {
    all_region_len = all_region_len + region[2] - region[1] + 1
  }
  cnv = all_reads_len/(all_region_len/50000)
  return(cnv)

}
