#' Computes 'proportions of tumor-like sites' PTS for each sample and creates a summary table
#'
#' @param dms_table A DataFrame specifying the chromosomal and genomic positions of selected DMS. It must include at least four columns: `dms_id` (a unique identifier for each DMS), `chr` (chromosome, e.g., 'chr1', 'chr2', etc.), `pos` (genomic position), and `type` ('hypo' for hypomethylated CpGs or 'hyper' for hypermethylated CpGs).
#' @param beta_table A DataFrame created using the `create_dms_beta_table` function, containing beta values (proportions) for the selected DMS (row names) for each sample (column names).Beta-values are calculated considering only reads with an alpha value of 1 that cover at least the minimum required number of CpGs (see `create_dms_beta_table` function).
#'
#' @return A DataFrame ('PTS Table') containing 'proportions of tumor-like sites' (PTS values) for each sample. It includes `PTS_hyper` (for hypermethylated DMS), `PTS_hypo` (for hypomethylated DMS), and `PTS_all` (combining both hyper- and hypomethylated DMS). The `PTS_all` metric can be used as a proxy for sample's tumor content (TC).
#' @export
#'
meter_quant_PTS <- function(dms_table, beta_table){

  assertthat::assert_that(all(c("dms_id", "chr", "pos", "type") %in% colnames(dms_table)),
                          msg = "The dms_table does not include required columns")

  dms_table <- as.data.frame(dms_table)
  beta_table <- as.data.frame(beta_table)


  ### select DMS of interest
  beta_table <- beta_table[which(rownames(beta_table) %in% dms_table$dms_id), ]


  ### hyper and hypo DMS of interest
  hyper_dms <- dms_table$dms_id[which(dms_table$type=='hyper')]
  hypo_dms <- dms_table$dms_id[which(dms_table$type=='hypo')]


  ### compute PTS of hyper DMS
  pts_hyper <- beta_table[which(rownames(beta_table) %in% hyper_dms), ]

  pts_hyper <- apply(pts_hyper, 2, function(i) {
    sum(i==100, na.rm = T)/(sum(i==100, na.rm = T) + sum(i==0, na.rm = T))
  })


  ### compute PTS of hypo DMS
  pts_hypo <- beta_table[which(rownames(beta_table) %in% hypo_dms), ]

  pts_hypo <- apply(pts_hypo, 2, function(i) {
    sum(i==0, na.rm = T)/(sum(i==100, na.rm = T) + sum(i==0, na.rm = T))
  })


  ### compute PTS of hyper or hypo DMS
  tmp_hypo <- beta_table[which(rownames(beta_table) %in% hypo_dms), ]
  tmp_hypo <- apply(tmp_hypo, 2, function(i) {
    sum(i==0, na.rm = T)
  })

  tmp_hyper <- beta_table[which(rownames(beta_table) %in% hyper_dms), ]
  tmp_hyper <- apply(tmp_hyper, 2, function(i) {
    sum(i==100, na.rm = T)
  })

  tmp_all <- apply(beta_table, 2, function(i) {
    (sum(i==100, na.rm = T) + sum(i==0, na.rm = T))
  })

  pts_all <- (tmp_hyper + tmp_hypo)/tmp_all

  pts_tab <- as.data.frame(cbind(pts_hyper, pts_hypo, pts_all))

  pts_tab$samp_id <- rownames(pts_tab)

  rownames(pts_tab) <- NULL

  pts_tab <- pts_tab[, c("samp_id", "pts_hyper", "pts_hypo", "pts_all")]

  return(pts_tab)


}



