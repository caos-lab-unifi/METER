
filter_cov_alpha100 <- function(path_bismark2bedGraph, path_cpg_file, path_read_table, path_out = NULL, min_sites = 6, remove_cpg=FALSE){

  if (is.null(path_out)) {
    path_out <- file.path(dirname(path = path_cpg_file), 'filter_cov_alpha100')
    dir.create(path = path_out, showWarnings = F)
  }

  read_table <- readRDS(path_read_table)
  cpg_file <- data.table::fread(path_cpg_file, data.table = F)


  ### check input
  read_table <- as.data.frame(read_table)
  cpg_file <- as.data.frame(cpg_file)


  ### FILTER READS
  tokeep <- read_table$seq_id[which(read_table$n_sites>=min_sites &
                                   (read_table$meth_perc==1 | read_table$meth_perc==0))]

  rm(read_table)

  ### FILTER cpg file
  cn_1 <- colnames(cpg_file)[1]

  ## create tmp colnames
  colnames(cpg_file) <- paste0('tmp_', 1:ncol(cpg_file))

  cpg_file <- cpg_file[which(cpg_file$tmp_1 %in% tokeep), ]
  rownames(cpg_file) <- NULL

  ## restore original colnames
  colnames(cpg_file) = c(cn_1, '', '', '', '')

  ### save filtered cpg file
  data.table::fwrite(cpg_file,
                     file = file.path(path_out, paste0('CpG_context_', gsub('.rds', '.txt.gz', basename(path_read_table), fixed = T))),
                     sep = '\t', compress = 'auto'
                     )


  ### run bismark bismark2bedGraph to create new cov file from the filtered CpG file
  in_file=file.path(path_out, paste0('CpG_context_', gsub('.rds', '.txt.gz', basename(path_read_table), fixed = T)))
  id=gsub('.rds', '', basename(path_read_table))

  system(paste(path_bismark2bedGraph, in_file, '--dir', path_out, '--output', id))

  if (remove_cpg==TRUE) {
    system(paste0('find ', path_out, ' -mindepth 1 -maxdepth 1 -type f -name "CpG_context*\\.txt\\.gz" -delete'))
    }

}



