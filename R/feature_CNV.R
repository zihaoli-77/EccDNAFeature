#' @title Get CNV
#' @description
#' A function, is used to get CNV.
#'
#' @param bamfile bamfile which you want to read
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
feature_CNV <- function(bamfile){

  gene_name <- unlist(strsplit(bamfile, "\\."))[1]

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
  bam_data <- readGAlignments(bamfile, use.names=TRUE, param=param1)

  # 初始化变量
  regions <- list()
  current_region <- NULL
  all_reads_len = 0
  # 遍历 BAM 数据
  for (i in seq_along(bam_data)) {
    start_pos <- start(bam_data[i])
    end_pos <- end(bam_data[i])
    all_reads_len = all_reads_len + bam_data[i]@elementMetadata@listData$AS

    if (is.null(current_region)) {
      # 初始化第一个区域
      current_region <- c(start_pos, end_pos)
    } else if (start_pos <= current_region[2]) {
      # 如果有重叠，更新当前区域的结束位置
      current_region[2] <- max(current_region[2], end_pos)
    } else {
      # 如果没有重叠，保存当前区域并开始新的区域
      regions <- append(regions, list(current_region))
      current_region <- c(start_pos, end_pos)
    }
  }

  # 添加最后一个区域
  if (!is.null(current_region)) {
    regions <- append(regions, list(current_region))
  }

  all_region_len = 0
  for (region in regions) {
    all_region_len = all_region_len + region[2] - region[1] + 1
  }
  cnv = all_reads_len/(all_region_len/50000)
  data <- data.frame(gene_name,cnv)
  return(data)

}
