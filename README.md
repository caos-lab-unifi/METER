Getting started with METER workflow
================
Marta Paoli
2025-01-24

# METER workflow tutorial

METER (METhylome analyzER) is a computational tool to analyze the
circulating tumor DNA (ctDNA) exploiting differentially methylated sites
(DMS) and regions (DMR) from low-pass whole genome bisulphite sequencing
(lpWGBS) data of cell-free DNA (cfDNA) samples. It comprises three
modules:

- **METER-quant** to measure the ctDNA, based on tumor-specific DMS

- **METER-detect** to classify cfDNA samples as ctDNA+ (METER+) or
  ctDNA- (METER-) (that is ctDNA is detected or not), based on
  tumor-specific DMR

- **METER-subtype** to infer specific subtype from ctDNA, based on tumor
  subtype-specific DMR

Each module relies on task-specific informative DMR and DMS identified
through the differential analysis between DNA-methylation data from
tumor-type specific tissue samples and control samples such as cell-free
DNA (cfDNA) samples from healthy donors or whole blood.

# METER-quant

In this module, the “proportion of tumor-like sites” (PTS) by sample is
measured from lpWGBS as a proxy of sample’s TC. PTS for each sample is
computed as the ratio of **fully methylated (beta=100%)
hypermethylated** or **fully unmethylated (beta=0%) hypomethylated**
**DMS** (that is DMS supporting a tumor-like methylation signal) to the
total fully methylated/unmethylated DMS covered. To minimize the effect
of bisulfite conversion errors and increase signal specificity, only
reads with alpha value of 100%, that is reads showing only methylated or
unmethylated CpG sites, and covering a minimum number of CpG sites (6
CpG sites by default) are considered for the computation of the beta
values.

## Functions:

The METER package provides a set of functions that use **Bismark (1)
“CpG_context”** files as a starting point to generate the Beta Table
necessary for computing PTS.

### create_read_table

Starting from Bismark (1) “CpG_context” file for a specific sample, this
function generates a table that summarizes information for each
sequencing read.

- **Input data and parameters**:
  - cpg_file: Bismark (1) “CpG_context” file for the sample analyzed,
    obtained using *–comprehensive* option (as strand-specific
    methylation is not of interest, see Bismark (1) manual
    <https://felixkrueger.github.io/Bismark/bismark/>)
  - dmr_table (optional, required for METER-detect module): a DataFrame
    reporting the chromosomal and genomic positions of selected DMR,
    containing at least 4 columns: “dmr_id” (unique identifier for a
    specific DMR), “chr” (chromosome reported as “chr1”, “chr2”…),
    “start” (starting genomic position), “end” (ending genomic
    position), and “type” (“hypo” for hypomethylated DMR or “hyper” for
    hypermethylated DMR). A data frame listing breast cancer-specific
    DMR for tumor content detection is included within METER package
    (dmr_filtered.rds).
- **Output**:
  - Read Table: a DataFrame containing summarized information for each
    sequencing read, including the position of the first CpG site within
    the read (“min_cpg_pos”), the position of the last CpG within the
    read (“max_cpg_pos”), the number of methylated CpGs (“n_meth”), the
    number of unmethylated CpGs (“n_unmeth”), and the total number of
    CpGs covered (“n_sites”). Additionally, if a DMR Table is provided,
    it specifies, for each read, the overlapping DMR (“dmr_id”) and the
    type of the overlapping DMR (“hypo” or “hyper”), where “overlapping”
    means the read is entirely contained within an DMR.

Note 1: this step takes a while, so if one has a lot of samples to
process parallelization is recommended when possible

Note 2: the file name of the output Read Table has to be chosen by the
user when saving the table, and it will serve as sample ID for
subsequent steps

``` r
path_cpg_files='/path/to/Bismark_files/'
path_out='/path/to/output_folder/create_read_table_noDMR/'

dir.create(path_out, showWarnings = F, recursive = T)

cpg_files=list.files(path_cpg_files, full.names = T, pattern = 'CpG_context')

# i=cpg_files[1]
parallel::mclapply(cpg_files, function(i){

  samp_id=gsub('CpG_context_', '', strsplit(basename(i), "_bismark_", fixed = T)[[1]][1])

  read_table=create_read_table(path_cpg_file = i)

  saveRDS(read_table, file.path(path_out, paste0(samp_id, '.rds')))

  }, mc.cores = 6)

rm(list = ls())


### take a look to one of the tables created
example_table=readRDS('/path/to/output_folder/create_read_table_noDMR/samp1.rds')
head(example_table)
```

### filter_cov_alpha100

Starting from Bismark (1) “CpG_context” file and the Read Table obtained
in the previous step for a specific sample, this function creates a
filtered “Bismark-like CpG_context” file, considering only reads with
alpha value of 100% (based on the samples’ Read Tables created in the
previous step) and covering a minimum number of CpGs (6 by default).
Subsequently, it employes the *bismark2bedGraph* function from Bismark
(1) to produce the corresponding coverage file (“.cov”), displaying for
each CpG the count of methylated/unmethylated reads covering its
position.

- **Input data and parameters**:
  - path_bismark2bedGraph: absolute path to *bismark2bedGraph* function
  - path_cpg_file: absolute path to Bismark (1) “CpG_context” file for
    the sample analyzed, obtained using the *–comprehensive* option (as
    strand-specific methylation is not of interest, see Bismark (1)
    manual <https://felixkrueger.github.io/Bismark/bismark/>)
  - path_read_table: absolute path to Read Table previously generated
  - path_out (optional): absolute path to output folder
  - min_sites (default min_sites=6): integer indicating the minimum
    number of CpGs a read should contain in order to be considered in
    the filtered “CpG_context” file and consequently in the filtered
    “coverage file”
  - remove_cpg (default remove_cpg=FALSE): TRUE if one wants to delete
    intermediate “CpG_context” file once the final “coverage file” has
    been created
- **Output**:
  - “Coverage file” generated by *bismark2bedGraph* function from
    Bismark (1), constructed considering only reads with alpha value of
    100% and covering a number of CpG sites ≥ min_sites.

Note 1: this step takes a while, so if one has a lot of samples to
process parallelization is recommended when possible

Note 2: The file name of the output .cov file derives from the base name
of the input Read Table (generated in the preceding step) that has been
chosen by the user

``` r
path_read_tables='/path/to/output_folder/create_read_table_noDMR/'
path_cpg_files='/path/to/Bismark_files/'
path_bis='/path/to/bismark2bedGraph'
out_folder='/path/to/output_folder/filter_cov_alpha100/'


### CpG files (output of bismark methylation extractor using --comprehensive option)
cpg_files=list.files(path_cpg_files, full.names = T, pattern = '^CpG_context*')

### bismark reads summary tables
read_tables=list.files(path_read_tables, full.names = T, pattern = '.rds')

### run filter_cov_alpha100 function on all input samples
# i=read_tables[1]
parallel::mclapply(read_tables, function(i){

  id=gsub('.rds', '', basename(i), fixed = T)

  path_cpg=cpg_files[which(grepl(id, cpg_files))]

  filter_cov_alpha100(path_bismark2bedGraph = path_bis, path_cpg_file = path_cpg, path_read_table = i, path_out = out_folder, remove_cpg = F)

}, mc.cores = 1)
```

### create_dms_beta_table

This function computes the beta values (as proportions) of selected DMS
for each sample, and creates a DataFrame with samples as column names
and CpG sites as row names. The computation of beta values is performed
using the filtered “coverage files” generated in the previous step.

- **Input data and parameters**:
  - dms_table: DataFrame reporting the chromosomal and genomic positions
    of selected DMS, containing at least 4 columns: “dms_id” (unique
    identifier for a specific DMS), “chr” (chromosome reported as
    “chr1”, “chr2”…), “pos” (genomic position), and “type” (“hypo” for
    hypomethylated CpGs or “hyper” for hypermethylated CpGs). A data
    frame listing breast cancer-specific DMS for tumor content
    quantification is included within METER package (dms_filtered.rds).
  - path_cov_files: absolute path to the folder that contains all the
    filtered Bismark (1) “coverage files” created for each sample in the
    previous step.
  - id_pattern (optional): A string to be used as input for the R
    `strsplit()` function to extract sample names from the base names of
    the input “coverage files.” If not specified (default = NULL), the
    sample names will be directly obtained from the base names of the
    input “coverage files”.
- **Output**:
  - Beta Table: a DataFrame containing the beta values (as proportions)
    of DMS across samples, with samples as column names and DMS as row
    names. Each DMS is reported with its unique “dms_id”, matching the
    identifiers provided in the DMS table. Beta values are computed
    considering only reads with alpha value of 100% and covering a
    minimum number of CpGs, as explained above.

Note 1: Samples’ IDs (column names of the output DMS beta table) are
derived from the base names of the input .cov files (that in turn derive
from the base names of the read tables). Bismark (1) attaches some
extensions to the .cov files generated that can be removed through the
“id_pattern” parameter.

``` r
path_cov_files='/path/to/output_folder/filter_cov_alpha100/'
path_dmss='/path/to/dms_filtered.rds'
path_out='/path/to/output_folder/'

dmss=readRDS(path_dmss)
head(dmss)

dms_beta_tab=create_dms_beta_table(dms_table = dmss, path_cov_files = path_cov_files, id_pattern = '.gz.bismark.cov.gz')

head(dms_beta_tab)

saveRDS(object = dms_beta_tab, file = file.path(path_out, 'dms_beta.rds'))
```

### meter_quant_PTS

computes PTS for each sample and creates a summary table

- **Input data and parameters**:
  - dms_table: a DataFrame reporting the chromosomal and genomic
    positions of selected DMS, containing at least 4 columns: “dms_id”
    (unique identifier for a specific DMS), “chr” (chromosome reported
    as “chr1”, “chr2”…), “pos” (genomic position), and “type” (“hypo”
    for hypomethylated CpGs or “hyper” for hypermethylated CpGs). A data
    frame listing breast cancer-specific DMS for tumor content
    quantification is included within METER package (dms_filtered.rds).
  - beta_table: a DataFrame as the one created in the previous step,
    containing the beta values (as proportions) of CpG sites across
    samples, with samples as column names and CpG sites as row names.
    Each CpG site must be reported with its unique “dms_id”, matching
    the identifiers provided in the DMS table. Beta values should be
    computed considering only reads with alpha value of 100% and
    covering a minimum number of CpGs, as explained above.
- **Output**:
  - PTS Table: a DataFrame that includes PTS values for hypermethylated
    DMS (PTS_hyper), hypomethylated DMS (PTS_hypo), and both hyper and
    hypomethylated DMS combined (PTS_all) for each sample. PTS_all
    measure can be used as a proxy of TC.

``` r
path_dmss='/path/to/dms_filtered.rds'
path_dms_beta='/path/to/output_folder/dms_beta.rds'
path_out='/path/to/output_folder/'


dmss=readRDS(path_dmss)
dms_beta=readRDS(path_dms_beta)

head(dmss)
head(dms_beta)

pts_tab=meter_quant_PTS(dms_table = dmss, beta_table = dms_beta)

head(pts_tab)

saveRDS(pts_tab, file.path(path_out, 'meter_quant_PTS.rds'))
```

# METER-detect

In this module, the “proportion of tumor-like sequenced reads” (PTR) is
computed from lpWGBS data for each sample, and a Z-Score approach, based
on the distribution of this measure in control samples (reference
model), is then used to classify samples as either ctDNA+ (METER+) or
ctDNA- (METER-). As in METER-quant module, to minimize the effect of
bisulfite conversion errors and increase signal specificity, only reads
with alpha value of 100% and covering a minimum number of CpGs (6 CpGs
by default) are considered for PTR computation. Specifically, PTR for
each sample is computed as the proportion of **fully methylated** reads
within the selected **hyper-DMR** and fully **unmethylated reads**
within the selected **hypo-DMR** (that is reads supporting a tumor-like
methylation signal) over the total reads with alpha=100% within the
selected DMR.

## Functions:

The METER package provides a set of functions that use **Bismark (1)
“CpG_context”** files as a starting point to generate the Read Tables
necessary for computing PTR.

### create_read_table

As outlined in the METER-quant section, this function utilizes a Bismark
(1) “CpG_context” file for a specific sample to generate a comprehensive
table summarizing information for each sequencing read. In detail, it
creates a DataFrame containing for each read the position of the first
CpG site within the read (“min_cpg_pos”), the position of the last CpG
within the read (“max_cpg_pos”), the number of methylated CpGs
(“n_meth”), the number of unmethylated CpGs (“n_unmeth”), and the total
number of CpGs covered (“n_sites”). Importantly, this module requires
the DMR Table to be provided. This allows the function to evaluate, for
each read, the overlapping DMR (“dmr_id”) and the type of the
overlapping DMR (“hypo” or “hyper”), where **“overlapping” means the
first and the last CpG sites within the read are both contained within a
DMR**. (See METER-quant module for details)

``` r
path_cpg_files='/path/to/Bismark_files/'
path_dmrs='/path/to/dmr_filtered.rds'
path_out='/path/to/output_folder/create_read_table/'


cpg_files=list.files(path_cpg_files, full.names = T, pattern = 'CpG_context')

dmrs=readRDS(path_dmrs)
head(dmrs)

# i=cpg_files[1]
parallel::mclapply(cpg_files, function(i){

  samp_id=gsub('CpG_context_', '', strsplit(basename(i), "_bismark_", fixed = T)[[1]][1])

  read_table=create_read_table(path_cpg_file = i, dmr_table = dmrs)

  saveRDS(read_table, file.path(path_out, paste0(samp_id, '.rds')))

  }, mc.cores = 6)


### take a look to one of the tables created
example_table=readRDS('/path/to/output_folder/create_read_table/samp1.rds')
head(example_table)
```

### meter_detect_PTR

computes PTR for each sample and creates a summary table

- **Input data and parameters**:
  - path_read_tables: absolute path to the folder containing all the
    Read Tables created for each sample using *create_read_table*
    function (see previous step and METER-quant section).
  - min_sites (default min_sites=6): integer indicating the minimum
    number of CpGs a read should contain in order to be considered for
    PTR computation
  - ncores (default ncores=1): integer value representing the number of
    processor cores to use for parallel processing of samples.
- **Output**:
  - PTR Table: a DataFrame containing PTR values for hypermethylated DMR
    (PTR_hyper), hypomethylated DMR (PTR_hypo), and both hyper and
    hypomethylated DMR combined (PTR_all) for each sample. **PTR_all**
    measure, integrating information from both hyper and hypomethylated
    DMR, should be used to classify samples. Specifically, using this
    measure, each sample can be classified as ctDNA+ (METER+) if its
    Z-Score exceeds a certain threshold (e.g. Z-Score=3 has been used in
    the METER study (2)). For breast cancer samples, values
    corresponding to various Z-Scores thresholds based on the
    distribution of PTR in control samples (reference model) are
    provided within METER package (PTR_z_scores.rds).

Note 1: Sample ids in the “samp_id” column of the PTR Tables derive from
the base name of the input Read Table files (generated in the preceding
step) that have been chosen by the user

``` r
path_read_tabs='/path/to/output_folder/create_read_table/'
path_out='/path/to/output_folder/'


ptr_tab=meter_detect_PTR(path_read_tables = path_read_tabs, min_sites = 6, ncores = 4)

head(ptr_tab)

saveRDS(ptr_tab, file.path(path_out, 'meter_detect_PTR.rds'))
```

After generating the PTR Table for lpWGBS samples, **PTR_all** measure
can be used to categorize samples as either ctDNA-positive or
ctDNA-negative. For breast cancer samples, we can use the
“PTR_z_score.rds” table available in the METER package. Specifically,
the values associated with PTR_all in the “est_method” column serve as
the benchmark thresholds. In the following example, a threshold
corresponding to a Z-Sore of 3 is employed.

``` r
path_zscore_table='/path/to/PTR_z_scores.rds'

zscore_table=readRDS(path_zscore_table)

zscore_table

thr=zscore_table$z3[which(zscore_table$est_method=='PTR_all')]

### classify samples
ptr_tab$meter_detected=ifelse(ptr_tab$ptr_all>=thr, TRUE, FALSE)

head(ptr_tab)
```

# METER-subtype

This module employes the Robust Partial Correlation (RPC) deconvolution
method, as implemented in the EpiDISH (3) R package (PMID: 28193155), to
evaluate tumor subtype in circulation. For this method, a reference
model is necessary, which is derived by computing the median beta value
for each selected DMR across reference samples representing the
components to be deconvoluted. Naturally, it is crucial that the DMR
included in the reference model are tumor subtype-specific, meaning they
should exhibit differential methylation between tumor subtype 1 and
tumor subtype 2, as well as between tumor and healthy cfDNA (details for
selecting such DMR are provided in the “Methods” of the METER study
(2)). Additionally, the RPC method requires a beta table containing beta
values for each selected DMR across all samples being analyzed. For
breast cancer samples, the METER package provides a reference model
including ER+/Her2-, TNBC, and healthy cfDNA components for 1605
subtype-specific DMR (dmr_filtered_subtype.rds).

## Functions:

### create_dmr_beta_table

This function computes beta by DMR values (as proportions) of selected
DMR for each sample, and creates a DataFrame with samples as column
names and DMR as row names. Beta values for each DMR are computed by
treating all CpGs within a specific DMR collectively as a single CpG.
Therefore, the beta value for each DMR is calculated based on the ratio
of methylated signals to the total methylation signals for all CpGs
within that DMR. In other words, it reflects the proportion of times all
CpGs within the DMR are observed as methylated by sequence reads.

- **Input data and parameters**:
  - dmr_table: a DataFrame reporting the chromosomal and genomic
    positions of selected subtype-specific DMR, containing at least 4
    columns: “dmr_id” (unique identifier for a specific DMR), “chr”
    (chromosome reported as “chr1”), “start” (starting genomic
    position), “end” (ending genomic position), and “type” (“hypo” for
    hypomethylated DMR or “hyper” for hypermethylated DMR). For breast
    cancer samples, the METER package provides a DataFrame including
    ER+Her2-, TNBC, and healthy cfDNA components for 1605
    subtype-specific DMR (dmr_filtered_subtype.rds).
  - path_cov_files: absolute path to the folder that contains all
    Bismark (1) “coverage files” for each sample.
  - id_pattern: A string to be used as input for the R `strsplit()`
    function to extract sample names from the base names of the input
    “coverage files.” If not specified (default = NULL), the sample
    names will be directly obtained from the base names of the input
    “coverage files”.
  - min_sites (default min_sites=0): minimum number of CpG sites that a
    DMR must include within a single sample. In simpler terms, it sets
    the smallest count of CpGs that need to be part of a DMR for it to
    be considered valid in the analysis of that sample.
- **Output**:
  - Beta Table: a DataFrame containing the beta by DMR values (as
    proportions) of selected subtype-specific DMR across samples, with
    samples as column names and DMR as row names. Each DMR is reported
    with its unique “dmr_id”, matching the identifiers provided in the
    DMR table.

Note 1: in this step it is not necessary to consider reads with alpha
value of 100% only to compute beta by DMR, therefore the following
example utilizes “original” Bismark (1) “coverage files” containing
information from all sequenced reads.

Note 2: in this deconvolution step it is recommended to consider the
beta of a specific DMR for a specific sample when the specific DMR for
the specific sample covers a minimum number of CpGs, so that the beta
computed is reliable (in our example we consider DMR containing a
minimum of 10 CpGs per sample)

``` r
path_cov_files='/path/to/Bismark_files/'
path_dmrs='/path/to/dmr_filtered_subtype.rds'
path_out='/path/to/output_folder/'

dmrs=readRDS(path_dmrs)
head(dmrs)

dmr_beta=create_dmr_beta_table(dmr_table = dmrs, path_cov_files = path_cov_files, id_pattern = '_bismark_', min_sites = 10)
head(dmr_beta)

saveRDS(object = dmr_beta, file = file.path(path_out, 'dmr_beta_subtype.rds'))
```

### meter_subtype

This function utilizes the Robust Partial Correlation (RPC)
deconvolution method, as implemented in the EpiDISH (3) R package (PMID:
28193155), to infer the tumor subtype from lpWGBS data of cfDNA samples.
The method estimates the proportions of tumor subtypes predefined in a
reference mwethylation matrix, which contains the median beta values, of
selected subtype-specific DMR, computed across subtype reference
samples. The reference samples included in the matrix represent the
components to be deconvoluted. Importantly, the components estimated
using EpiDISH (3) are designed to sum to 1, as the input reference
matrix is expected to account for all possible components of a sample.
Therefore, the reference matrix must include a column labeled
“healthy_cfDNA”, representing the healthy cell-free DNA (cfDNA)
component, alongside columns representing the tumor subtypes of
interest, as METER-subtype is specifically designed for cfDNA samples.

- **Input data and parameters**:
  - dmr_beta: a DataFrame (as obtained through the create_dmr_beta_table
    function in the previous step) containing the beta by DMR values (as
    proportions) of selected subtype-specific DMR across samples to be
    inspected, with samples as column names and DMR as row names. Each
    DMR must be reported with a unique “dmr_id”.
  - ref_mat: a matrix or DataFrame containing the median beta values
    computed across subtype reference samples for each subtype-specific
    DMR, with subtypes as column names and DMR as row names. Each DMR
    must be reported with its unique “dmr_id”, matching the identifiers
    provided in the dmr_beta table. The column names of the matrix
    represent the components to be deconvoluted, and must include a
    column labeled “healthy_cfDNA” representing the healthy cfDNA
    component.
- **Output**:
  - a DataFrame containing for each sample (row names) the estimated
    proportion of each component included in ref_mat (column names). An
    additional column labeled “pred_subtype”, alongsides columns
    relative to each component, reports the inferred subtype based on
    the predominant tumor subtype (not healthy cfDNA) estimated
    component.

``` r
path_subtype_ref='/path/to/ERpos_TNBC_ref.rds'
path_dmr_beta_sub='/path/to/output_folder/dmr_beta_subtype.rds'
path_out='/path/to/output_folder/'

dir.create(path_out, showWarnings = F)

subtype_ref=readRDS(path_subtype_ref)
dmr_beta_sub=readRDS(path_dmr_beta_sub)

est_sub=meter_subtype(ref_mat = subtype_ref, dmr_beta = dmr_beta_sub)

saveRDS(est_sub, file.path(path_out, 'RPC_components.rds'))
```

# References

- 1)  Krueger F, Andrews SR. Bismark: a flexible aligner and methylation
      caller for Bisulfite-Seq applications. Bioinformatics. 2011 Jun
      1;27(11):1571-2. doi: 10.1093/bioinformatics/btr167. Epub 2011
      Apr 14. PMID: 21493656; PMCID: PMC3102221.

- 2)  doi: <https://doi.org/10.1101/2024.06.10.598204>

- 3)  Teschendorff AE, Breeze CE, Zheng SC, Beck S. A comparison of
      reference-based algorithms for correcting cell-type heterogeneity
      in Epigenome-Wide Association Studies. BMC Bioinformatics. 2017
      Feb 13;18(1):105. doi: 10.1186/s12859-017-1511-5. PMID: 28193155;
      PMCID: PMC5307731.
