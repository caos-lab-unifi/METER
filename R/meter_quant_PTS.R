
meter_quant_PTS <- function(dms_table, beta_table){

  dms_table <- as.data.frame(dms_table)
  beta_table <- as.data.frame(beta_table)


  ### select DMS of interst
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



