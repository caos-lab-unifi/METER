#' Compute PTR of a sample considering only reads with alpha=100
#' Create PTR table for all study samples
#'
#'
#' input:
#' - dmr table with specific columns: dmr_id, type (hyper/hypo)
#' - reads-dmr summary table with specific columns: read_id, dmr_id, meth_perc (alpha of read, i.e. degree of methylation of read)
#'
#' NB: dmr_id in the reads-dmr table sholud be a subset of the dmr_id in the dmr table
#'
#' output: data.frame with 1 row containing sample_id, ptr hyper, ptr hypo and ptr glob
#'
#'
#'

create_ptr_table <- function(path_read_tables, min_sites = 6, ncores = 1){

  ### check input
  system_cores <- parallel::detectCores()

  if (!is.na(system_cores)){
    assertthat::assert_that(ncores < system_cores)
  }

  ### list files
  lf=list.files(path_read_tables, full.names = T, pattern = '.rds')


  ### compute PTR by sample
  ptr_data <- parallel::mclapply(mc.cores = ncores, lf, function(i) {

    read_dmr_table <- readRDS(i)

    id=gsub('\\.rds', '', basename(i))


    ### check input
    read_dmr_table <- as.data.frame(read_dmr_table)
    # assertthat::assert_that(is.character(id))
    # assertthat::assert_that(c('dmr_id') %in% colnames(dmr_table) & is.character(dmr_table$dmr_id))
    # assertthat::assert_that(c('chr') %in% colnames(dmr_table))
    # assertthat::assert_that(c('start') %in% colnames(dmr_table) & is.integer(dmr_table$start))
    # assertthat::assert_that(c('end') %in% colnames(dmr_table) & is.integer(dmr_table$end))
    # assertthat::assert_that(nrow(read_dmr_table)>0)
    # assertthat::assert_that(nrow(dmr_table)>0)


    ### prepare read_dmr_table
    ## select columns
    read_dmr_table <- read_dmr_table[, c("seq_id", "n_sites", "meth_perc", "dmr_id", "dmr_type")]

    ## remove rows with NA, i.e. reads not overlapping iDMR
    read_dmr_table <- read_dmr_table[complete.cases(read_dmr_table), ]

    ## select reads with alpha=100
    read_dmr_table <- read_dmr_table[which(read_dmr_table$meth_perc==1 |
                                             read_dmr_table$meth_perc==0), ]

    rownames(read_dmr_table) <- NULL


    ### check if reads are left
    # assertthat::assert_that(nrow(read_dmr_table)>0)


    ### compute proportion of tumor-like reads (PTR)
    read_dmr_table$tumor_like <- dplyr::case_when(
      read_dmr_table$meth_perc==1 & read_dmr_table$dmr_type=='hyper' ~ 1,
      read_dmr_table$meth_perc==0 & read_dmr_table$dmr_type=='hypo' ~ 1,
      !is.na(read_dmr_table$meth_perc) & !is.na(read_dmr_table$dmr_type) ~ 0
    )

    ptr_hyper=sum(read_dmr_table$tumor_like==1 & read_dmr_table$dmr_type=='hyper' & read_dmr_table$n_sites>=min_sites) /
      sum(read_dmr_table$dmr_type=='hyper' & read_dmr_table$n_sites>=min_sites)

    ptr_hypo=sum(read_dmr_table$tumor_like==1 & read_dmr_table$dmr_type=='hypo' & read_dmr_table$n_sites>=min_sites) /
      sum(read_dmr_table$dmr_type=='hypo' & read_dmr_table$n_sites>=min_sites)

    ptr_all=sum(read_dmr_table$tumor_like==1 & read_dmr_table$n_sites>=min_sites) /
      sum(read_dmr_table$n_sites>=min_sites)

    ptr_samp=data.frame(samp_id=id, ptr_hyper=ptr_hyper, ptr_hypo=ptr_hypo, ptr_all=ptr_all)

    return(ptr_samp)

  })


  ### create PTR table
  ptr_table=as.data.frame(dplyr::bind_rows(ptr_data))

  return(ptr_table)

}










