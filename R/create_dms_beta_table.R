#' Create iDMS beta table for study samples
#'
#'
#' input:
#' for each sample:
#' - iDMS of interest
#' - Bismark .cov files for all study samples created in during processing, created considering reads with alpha=100 only
#'
#' output:
#' iDMS beta table for study samples
#'
#'
#'



create_dms_beta_table <- function(dms_table, path_cov_files, id_pattern=NULL){


  ### list files
  lf <- list.files(path_cov_files, full.names = T, pattern = '.cov')


  ### DMS TABLE Granges
  gr_dms <- GenomicRanges::GRanges(seqnames = dms_table$chr,
                                ranges = IRanges::IRanges(start = dms_table$pos,
                                                          width = 1),
                                dms_id=dms_table$dms_id,
                                type=dms_table$type)


  ### create iDMS table
  beta_data=lapply(lf, function(i) {

    if (is.null(id_pattern)) {
      id <- basename(i)
    } else {
      id <- strsplit(basename(i), id_pattern, fixed = T)[[1]][1]
    }

    cov <- data.table::fread(i, header = F, data.table = F)
    colnames(cov) <- c('chr', 'start', 'end', 'beta', 'nC', 'nT')


    ### overlap iDMS with cov
    cov <- GenomicRanges::GRanges(seqnames = cov$chr,
                                  ranges = IRanges::IRanges(start = cov$start, width = 1),
                                  beta=cov$beta)

    cov <- IRanges::mergeByOverlaps(query = cov, subject = gr_dms, type='within')
    cov <- as.data.frame(cov)


    ### prepare final data
    cov <- cov[, c("gr_dms.dms_id", "cov.beta")]
    colnames(cov) <- c("dms_id", id)


    ### convert to datatable and setkey
    data.table::setDT(cov)
    data.table::setkey(cov)

    cov

  })

  beta_data <- as.data.frame(Reduce(function(x,y) {merge(x,y, by=c('dms_id'), all=T)}, beta_data))

  rownames(beta_data) <- beta_data$dms_id

  beta_data <- subset(beta_data, select = -c(dms_id))

  return(beta_data)

}



