#' @import data.table
create_read_table <- function(path_cpg_file, dmr_table=NULL){

  ### load Bismark "CpG_context" file
  cpg_file=fread(path_cpg_file, header = F)


  ### summarise Bismark CpG file per read: obtain dataframe with unique read per row
  colnames(cpg_file) <- c('seq_id', 'meth_state', 'chr', 'start', 'meth_call')

  setDT(cpg_file)

  cpg_file <- cpg_file[, .(
    min_cpg_pos = min(start),
    max_cpg_pos = max(start),
    n_meth = sum(meth_call == 'Z'),
    n_sites = .N), by = .(seq_id, chr)]

  cpg_file[, meth_perc := n_meth / (n_sites)]
  cpg_file[, n_unmeth := (n_sites - n_meth)]

  cpg_file=as.data.frame(cpg_file)


  if (is.null(dmr_table)) {

    return(cpg_file)

  } else {


    dmr_table <- as.data.frame(dmr_table)

    ### create dmr GRanges
    gr_dmr <- GenomicRanges::GRanges(seqnames = dmr_table$chr,
                                     ranges = IRanges::IRanges(start = dmr_table$start,
                                                               end = dmr_table$end),
                                     dmr_id=dmr_table$dmr_id,
                                     dmr_type=dmr_table$type)


    ### create reads GRanges
    gr_reads <- GenomicRanges::GRanges(seqnames = cpg_file$chr,
                                    ranges = IRanges::IRanges(start = cpg_file$min_cpg_pos,
                                                              end = cpg_file$max_cpg_pos),
                                    seq_id = cpg_file$seq_id,
                                    n_meth = cpg_file$n_meth,
                                    n_unmeth = cpg_file$n_unmeth,
                                    n_sites=cpg_file$n_sites,
                                    meth_perc=cpg_file$meth_perc)


    ### create final table showing if a specific read overlaps (maps within) a specific DMR
    ## NB: overlap between read and DMR is performed considering the first and the last CpG site within the read
    overlap <- IRanges::mergeByOverlaps(query = gr_reads, subject = gr_dmr, type='within') ## first and last CpG sites must be both within the DMR

    overlap <- as.data.frame(overlap)

    overlap <- overlap[, c("gr_reads.seq_id", "gr_reads.seqnames", "gr_reads.start", "gr_reads.end",
                           "gr_reads.n_meth", "gr_reads.n_unmeth", "gr_reads.n_sites", "gr_reads.meth_perc",
                           "gr_dmr.dmr_id", "gr_dmr.dmr_type")]

    colnames(overlap) <- gsub('gr_dmr.', '', colnames(overlap), fixed = T)
    colnames(overlap) <- gsub('gr_reads.', '', colnames(overlap), fixed = T)
    colnames(overlap)[which(colnames(overlap) == 'seqnames')] <- 'chr'
    colnames(overlap)[which(colnames(overlap) == 'start')] <- 'min_cpg_pos'
    colnames(overlap)[which(colnames(overlap) == 'end')] <- 'max_cpg_pos'

    ## reads NOT overlapping DMR (we want to keep these reads)
    out_reads <- cpg_file[which(!(cpg_file$seq_id %in% overlap$seq_id)), ]
    out_reads$dmr_id <- NA
    out_reads$dmr_type <- NA

    ## final table
    overlap <- rbind(overlap, out_reads)

    return(overlap)

  }

}





