# Scripts for SAMOSA-Tag analyses

## Scripts
Please see [`SAMOSA-Tag_processing.ipynb`](../notebooks/SAMOSA-Tag_processing.ipynb) for examples of how these scripts are used.

Scripts assume the following top-level directory `${TOP_DIR}/analyses/OS152/` with external files at paths noted below.

* [`01_compute_autocors_persample.py`](./01_compute_autocors_persample.py): Compute per-molecule accessibility autocorrelation profiles. These autocorrelograms are used to define the type of nucleosome positioning observed on each read.

* [`02_cluster_autocorrelograms_tn5.py`](./02_cluster_autocorrelograms_tn5.py): Cluster the single-molecule autocorrelation profiles produced by `01_compute_autocors_persample.py` to define _fiber types_ - consensus profiles of nucleosome positioning.

* [`03_samosa_feature_signal.py`](./03_samosa_feature_signal.py): Examine the enrichment of single-molecule accessibility at specific genomic sites. Here, we examine CTCF sites defined in U2OS cells (see `external datasets` below for more information)

* [`05_cpg2pickle.py`](./05_cpg2pickle.py): Convert 5mC methylation predictions produced by [`primrose`](https://github.com/PacificBiosciences/primrose) and stored in BAM tags ML and MM per-read into arrays for downstream analyses. 

* [`06_primrose_samosa_ctcf_integration_OS.py`](./06_primrose_samosa_ctcf_integration_OS.py): Examine whether 5mC methylation at CpG sites per molecule are enriched at U2OS CTCF sites, similar to `03_samosa_feature_signal.py`

* [`07_compare_cpg_samosa_OS_data.py`](./07_compare_cpg_samosa_OS_data.py): For each molecule, link together 5mC predictions produced by [`primrose`](https://github.com/PacificBiosciences/primrose) and accessbility profiles produced by the [SAMOSA-ChAAT computational pipeline](https://github.com/RamaniLab/SAMOSA-ChAAT)

* [`08_fishers_methylation_OS_SMRT_tag.py`](./08_fishers_methylation_OS_SMRT_tag.py): Stratify single molecules by feautres of 5mC methylation at CpGs (CpG density, methylation probability), and test for an enrichment of fiber types defined by `02_cluster_autocorrelograms_tn5.py`. 

* [`08b_reps_samosa_fishers_OS.py`](./08b_reps_samosa_fishers_OS.py): Stratify single molecules by feautres of 5mC methylation at CpGs (CpG density, methylation probability) as in `08_fishers_methylation_OS_SMRT_tag.py`, and test for an enrichment of fiber types defined by `02_cluster_autocorrelograms_tn5.py`, further stratified by sample replicates. 

* [`09_tss2endmatrix_pb.py`](./09_tss2endmatrix_pb.py): Examine whether the ends of SAMOSA-Tag reads are enriched in transcription start sites (TSS) of genes. 


## External datasets
GRCh38 reference genome:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
```

Scripts assume this file is at path `${TOP_DIR}/ref/GRCh38/hg38.fa`

U2OS ChIP-seq data
Data can be obtained from the following GEO entries:
* [`GSE87831`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87831)

## Software versions

* [`primrose v1.3.0`](https://github.com/PacificBiosciences/primrose/releases/tag/v1.3.0)

* `python >= 3.8`

These scripts also rely on the following python packages:

* [`numpy v1.19.0`](https://github.com/numpy/numpy/releases/tag/v1.19.0)

* [`pandas v1.2.4`](https://github.com/pandas-dev/pandas/releases/tag/v1.2.4)

* [`pysam v0.19.1`](https://github.com/pysam-developers/pysam/releases/tag/v0.19.1)

* [`scipy v1.9.0`](https://github.com/scipy/scipy/releases/tag/v1.9.0)

* [`scanpy v1.8.1`](https://github.com/scverse/scanpy/releases/tag/1.8.1)

* [`numba v0.55.0`](https://github.com/numba/numba/releases/tag/0.55.0)



