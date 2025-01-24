#' Create table of DMS beta-values
#'
#' @param dms_table A DataFrame specifying the chromosomal and genomic positions of selected DMS, with at least 4 required columns: `dms_id` (unique identifier for each DMS), `chr` (chromosome, e.g., 'chr1', 'chr2'), `pos` (genomic position), and `type` (methylation status, either 'hypo' for hypomethylated or 'hyper' for hypermethylated DMS)
#' @param path_cov_files The absolute path to the folder containing the filtered coverage files of samples to be analyzed, which should have been created using the `METER::filter_cov_alpha100` function.
#' @param id_pattern A string to be used as input for the R `strsplit()` function to extract sample names from the base names of the input "coverage files." If not specified (default = NULL), the sample names will be directly obtained from the base names of the input "coverage files".
#'
#' @return A DataFrame (Beta Table) containing the beta-values (proportions) of the selected DMS, with DMS as row names (identified by their unique `dms_id`, consistent with the `dms_table`) and samples as column names. Beta-values are calculated based on reads with alpha value of 1 and that cover a minimum number of CpGs.
#' @export
#'
#'
create_dms_beta_table <- function(dms_table, path_cov_files, id_pattern=NULL){

  assertthat::assert_that(all(c("dms_id", "chr", "pos", "type") %in% colnames(dms_table)),
                          msg = "The dms_table does not include required columns")


  ### list files
  lf <- list.files(path_cov_files, full.names = T, pattern = '\\.cov')

  assertthat::assert_that(length(lf)>0,
                          msg = "the coverage file folder is empty")


  ### DMS TABLE Granges
  gr_dms <- GenomicRanges::GRanges(seqnames = dms_table$chr,
                                ranges = IRanges::IRanges(start = dms_table$pos,
                                                          width = 1),
                                dms_id=dms_table$dms_id,
                                type=dms_table$type)


  ### create DMS table
  beta_data=lapply(lf, function(i) {

    if (is.null(id_pattern)) {
      id <- basename(i)
    } else {
      id <- strsplit(basename(i), id_pattern, fixed = T)[[1]][1]
    }

    cov <- data.table::fread(i, header = F, data.table = F)

    assertthat::assert_that(length(colnames(cov))==6, msg = "The coverage file is not in the right format")

    colnames(cov) <- c('chr', 'start', 'end', 'beta', 'nC', 'nT')


    ### overlap DMS with cov
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



