% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_read_table.R
\name{create_read_table}
\alias{create_read_table}
\title{Create a summary table of a sample's sequencing reads}
\usage{
create_read_table(path_cpg_file, dmr_table = NULL)
}
\arguments{
\item{path_cpg_file}{The absolute path to the 'CpG_context' file generated by Bismark for the analyzed sample. This file must be produced using the `bismark_methylation_extractor` function with the -`-comprehensive` option, as strand-specific methylation information is not required. The file should contain five columns reporting in the following order: sequencing read ID, methylation status, chromosomal position, genomic position, and methylation call of the CpG covered by the sequencing read. For more details, consult the Bismark manual: https://felixkrueger.github.io/Bismark/bismark/.}

\item{dmr_table}{A DataFrame specifying the chromosomal and genomic positions of selected DMR, with at least 4 required columns: `dmr_id` (unique identifier for each DMR), `chr` (chromosome, e.g., 'chr1', 'chr2'), `start` (starting genomic position), `end` (ending genomic position), and `type` (methylation status, either 'hypo' for hypomethylated or 'hyper' for hypermethylated DMR).}
}
\value{
A DataFrame (Read Table) summarizing each sequencing read included in the 'CpG_context' file of a sample, including the columns: `seq_id` (the sequence read identifier), `min_cpg_pos` (position of the first CpG within the read), `max_cpg_pos` (position of the last CpG within the read), `n_meth` (number of methylated CpGs), `n_unmeth` (number of unmethylated CpGs), and `n_sites` (total CpGs covered). If a DMR Table is provided, it also includes `dmr_id` (unique identifier of the overlapping DMR) and `type` ('hypo' for hypomethylated or 'hyper' for hypermethylated DMR), where overlapping means the read is fully contained within the DMR.
}
\description{
Create a summary table of a sample's sequencing reads
}
