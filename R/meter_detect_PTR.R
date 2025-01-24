#' Computes 'proportions of tumor-like reads' PTR for each sample and creates a summary table
#'
#' @param path_read_tables The absolute path to the folder containing all the Read Tables generated for each sample using the create_read_table function.
#' @param min_sites An integer specifying the minimum number of CpGs a read must contain to be included in the PTR calculation (default: min_sites = 6).
#' @param ncores An integer specifying the number of processor cores to use for parallel processing of samples (default: ncores = 1).
#'
#' @return A DataFrame ('PTR Table') containing PTR values for each sample, including `PTR_hyper` (for hypermethylated DMR), `PTR_hypo` (for hypomethylated DMR), and `PTR_all` (combining both hyper- and hypomethylated DMR). The `PTR_all` metric, which integrates information from both hyper- and hypomethylated DMR, is recommended for classifying samples as ctDNA+/-.
#' @export
#'
meter_detect_PTR <- function(path_read_tables, min_sites = 6, ncores = 1){

  # ### check input
  # system_cores <- parallel::detectCores()
  #
  # if (!is.na(system_cores)){
  #   assertthat::assert_that(ncores < system_cores)
  # }

  ### list files
  lf=list.files(path_read_tables, full.names = T, pattern = '.rds')


  ### compute PTR by sample
  ptr_data <- parallel::mclapply(mc.cores = ncores, lf, function(i) {

    read_dmr_table <- readRDS(i)

    id=gsub('\\.rds', '', basename(i))

    read_dmr_table <- as.data.frame(read_dmr_table)


    ### prepare read_dmr_table
    ## select columns
    read_dmr_table <- read_dmr_table[, c("seq_id", "n_sites", "meth_perc", "dmr_id", "dmr_type")]

    ## remove rows with NA, i.e. reads not overlapping DMR
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










