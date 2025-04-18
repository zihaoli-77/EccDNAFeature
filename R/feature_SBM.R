#' @title Get SBM
#' @description
#' A function, is used to get SBM.
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
#' @return SBM
#' @export
#'
#'
feature_SBM <- function(bamfile){


  param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                           hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                           isFirstMateRead = NA, isSecondMateRead = NA,
                                           isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                           isDuplicate = NA, isSupplementaryAlignment = NA),
                        what = "seq")
  primary <- readGAlignments(bamfile, param = param,use.names = T)


  a<-mcols(primary)
  a.seqnames <- seqnames(primary)
  a.start<-start(primary)
  a.end<-end(primary)
  a.strand<-strand(primary)
  a.cigar <- cigar(primary)
  df1<-(a)
  a.id<-rownames(df1)
  df1$p.chr <- a.seqnames
  df1$p.start<-a.start
  df1$p.end<-a.end
  df1$p.strand<-a.strand
  df1$p.cigar <- a.cigar
  df1$p.signal <- NA
  df1$p.signal <- ifelse(grepl("^[0-9]+M.*[0-9]+[SH]$|^[0-9][SH][0-9]+M", df1$p.cigar), "A","B")
  df1$id<-a.id




  hg.38 <- BSgenome.Hsapiens.UCSC.hg38

  # 生成非标准染色体名称时强制转为字符型
  non_grch38_chr <- as.character(setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY")))

  # 过滤前转换 p.chr 列为字符型
  df1$p.chr <- as.character(df1$p.chr)
  df1 <- df1[!(df1$p.chr %in% non_grch38_chr), ]

  # non_grch38_chr <- setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY"))
  #
  # df1 <- df1[!(df1$p.chr %in% non_grch38_chr),]



  qname1<- df1$id
  start1<-df1$p.start
  end1 <- df1$p.end
  p.signal <- df1$p.signal
  ##
  pA.qname <- qname1[p.signal == 'A']
  #df1$p.chr <- as.factor(df1$p.chr)
  df1$p.chr <- as.character(df1$p.chr)
  #pA.chr <- df1$p.chr[p.signal == 'A']
  pA.chr <- as.character(df1$p.chr[p.signal == 'A'])

  pA.end <- end1[p.signal == 'A']
  pB.qname <- qname1[p.signal == 'B']
  #pB.chr <- df1$p.chr[p.signal == 'B']
  pB.chr <- as.character(df1$p.chr[p.signal == 'B'])
  pB.start <- start1[p.signal == 'B']
  gr1 <- GRanges(seqnames = pA.chr,ranges = IRanges(start = pA.end-1,end = pA.end+2))
  gr2 <- GRanges(seqnames = pB.chr,ranges = IRanges(start = pB.start-2,end = pB.start+1))
  pa.seq<-getSeq(hg.38,gr1)
  pb.seq<-getSeq(hg.38,gr2)
  pa.seq.df <- as.character(pa.seq)
  pb.seq.df <- as.character(pb.seq)

  pa.df<- as.data.frame(pA.qname)
  pa.df$chr.1 <- pA.chr
  pa.df$motif.1 <- pa.seq.df
  names(pa.df) <- c("Qname","Chr.1","Motif.1")

  pb.df <- as.data.frame(pB.qname)
  pb.df$chr.1 <- pB.chr
  pb.df$motif.1 <- pb.seq.df
  names(pb.df) <- c("Qname","Chr.1","Motif.1")

  one <- rbind(pa.df,pb.df)

  one <- one[!grepl("N", one$Motif.1), ]

  sbpm<- as.data.frame(prop.table(table(one$Motif.1)))
  names(sbpm) <- c('motif','frequency')
  number_motif<-"4"
  base_pairs <- c("A", "T", "C", "G")
  all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
  all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
  empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
  names(empty.df)<-'motif'
  sBPM<-merge(empty.df,sbpm,by='motif',all.x=TRUE)
  sBPM[is.na(sBPM)] <- 0

  return(sBPM)

}
