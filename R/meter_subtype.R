#' Infer tumor subtype from cfDNA samples
#'
#' @param ref_mat A matrix or DataFrame of median beta-values for DMR (row names) for reference components (column names), including a column labeled `healthy_cfDNA` with median beta-values of the reference cfDNA component.
#' @param dmr_beta A DataFrame containing beta-values for DMR (row names) across cfDNA samples to be analyzed (column names).
#'
#' @return A DataFrame where each row corresponds to an analyzed sample, containing the proportions of each component and an additional `pred_subtype` column indicating the inferred subtype.
#' @export
#'
meter_subtype <- function(ref_mat, dmr_beta){

  common_dmr=intersect(rownames(ref_mat), rownames(dmr_beta))

  assertthat::assert_that("healthy_cfDNA" %in% colnames(ref_mat), msg = "'healthy_cfDNA' column is missing in ref_mat")
  assertthat::assert_that(length(common_dmr)>0, msg = "There are no commom dmr ids between ref_mat and dmr_beta")

  ref_mat=ref_mat[common_dmr, ]

  dmr_beta=dmr_beta[common_dmr, ]

  ref_mat=ref_mat[order(row.names(ref_mat)), ]
  dmr_beta=dmr_beta[order(row.names(dmr_beta)), ]

  ref_mat=as.matrix(ref_mat)
  dmr_beta=as.matrix(dmr_beta)


  ### Run EpiDISH
  out_epidish=EpiDISH::epidish(beta.m = dmr_beta, ref.m = ref_mat, method = "RPC")$estF

  out_epidish=as.data.frame(out_epidish)

  sub_out=out_epidish[ , -which(colnames(out_epidish) == "healthy_cfDNA")]

  out_epidish$pred_subtype = apply(sub_out, 1, function(row) {
    if (all(row == row[1])) {
      NA_character_
    } else {
      colnames(sub_out)[which.max(row)]
    }
  })

  rm(sub_out)

  return(out_epidish)

}








