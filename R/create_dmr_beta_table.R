#' Create table of DMR beta-values
#'
#' @param dmr_table A DataFrame specifying the chromosomal and genomic positions of selected DMR, with at least 4 required columns: `dmr_id` (unique identifier for each DMR), `chr` (chromosome, e.g., 'chr1', 'chr2'), `start` (starting genomic position), `end` (ending genomic position), and `type` (methylation status, either 'hypo' for hypomethylated or 'hyper' for hypermethylated DMR).
#' @param path_cov_files The absolute path to the folder containing the coverage files for the samples to be analyzed, with a .cov extension. These files must follow the formatting of Bismark coverage files (see: https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
#' @param id_pattern A string to be used as input for the R `strsplit()` function to extract sample names from the base names of the input "coverage files." If not specified (default = NULL), the sample names will be directly obtained from the base names of the input "coverage files".
#' @param min_sites An integer indicating the minimum number of CpG sites that a DMR must include within a single sample (default = 0).
#'
#' @return A DataFrame (Beta Table) containing the beta-values (proportions) of the selected DMR, with DMR as row names (identified by their unique `dmr_id`, consistent with the `dmr_table`) and samples as column names.
#' @export
#'

create_dmr_beta_table <- function(dmr_table, path_cov_files, id_pattern=NULL, min_sites=0){

  ### list files
  lf <- list.files(path_cov_files, full.names = T, pattern = '.cov')


  ### DMR TABLE Granges
  gr_dmr <- GenomicRanges::GRanges(seqnames = dmr_table$chr,
                                ranges = IRanges::IRanges(start = dmr_table$start,
                                                          end = dmr_table$end),
                                dmr_id=dmr_table$dmr_id,
                                type=dmr_table$type)


  ### create DMR table
  beta_data <- lapply(lf, function(i) {

    if (is.null(id_pattern)) {
        id <- basename(i)
      } else {
        id <- strsplit(basename(i), split = id_pattern, fixed = T)[[1]][1]
      }

    cov <- data.table::fread(i, header = F, data.table = F)
    colnames(cov)<- c('chr', 'start', 'end', 'beta', 'nC', 'nT')


    ### overlap DMR with cov
    cov <- GenomicRanges::GRanges(seqnames = cov$chr,
                                  ranges = IRanges::IRanges(start=cov$start, width=1),
                                  nC = cov$nC,
                                  nT = cov$nT,
                                  beta=cov$beta)

    overlap <- IRanges::mergeByOverlaps(query = cov, subject = gr_dmr, type='within')
    overlap <- as.data.frame(overlap)


    ### compute beta by region and number of CpG sites by region
    overlap <- dplyr::group_by(.data = overlap, gr_dmr.dmr_id)
    overlap <- dplyr::summarise(.data = overlap,
                                nC=sum(cov.nC),
                                nT=sum(cov.nT),
                                nsites=dplyr::n())

    overlap$beta <- overlap$nC/(overlap$nC + overlap$nT)

    overlap <- as.data.frame(overlap)


    ### set beta by DMR to NA when DMR comprises less than min_sites
    overlap$beta <- dplyr::case_when(overlap$nsites>=min_sites ~ overlap$beta,
                                     overlap$nsites<min_sites ~ NA_real_)


    ### final data
    overlap <- overlap[, c("gr_dmr.dmr_id", "beta")]
    colnames(overlap) <- c('dmr_id', id)


    ### convert to datatable and setkey
    data.table::setDT(overlap)
    data.table::setkey(overlap)

    return(overlap)

  })

  beta_data=as.data.frame(Reduce(function(x,y) {merge(x,y, by=c('dmr_id'), all=T)}, beta_data))

  rownames(beta_data) <- beta_data$dmr_id

  beta_data <- subset(beta_data, select = -c(dmr_id))

  return(beta_data)

}








