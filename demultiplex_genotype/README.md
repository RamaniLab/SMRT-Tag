# HG002, HG003, and HG004 genotype demultiplexing

## Scripts
Please see [`SMRT-Tag_processing.ipynb`](../notebooks/SMRT-Tag_processing.ipynb) for examples of how these scripts are used.

Scripts assume the following top-level directory: `${TOP_DIR}/analyses/HG/demultiplex_genotype/` with external files at paths noted below.

* [`demux_process_truth_sets.sh`](./demux_process_truth_sets.sh): create calling regions BED and private sites VCF files.

* [`demux_plp2vcf.sh`](./demux_process_truth_sets.sh): Naive variant calling using `samtools mpileup` and [`plp2vcf.py`](./plp2vcf.py), which is assumed to be in `$PATH`. Input is a BAM file aligned with `pbmm2`/`minimap2` using same parameters as in [`variant_calling`](../variant_calling)). 

* [`demux_compare.sh`](./demux_compare.sh): Compare variant calls in benchmark regions.

## External datasets
GRCh37 hs37d5 reference genome:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh37/hs37d5.fa.gz
```
Scripts assume this file is at path `${TOP_DIR}/ref/GRCh37/hs37d5.fa`

GIAB GRCh37 v3.0 genome stratification files:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/v3.0-stratifications-GRCh37.tar.gz
```
Scripts require just one stratification file from the directory at the following path:
`${TOP_DIR}/ref/GRCh37/union/GRCh37_notinalldifficultregions.bed.gz`

GIAB v4.2.1 benchmark call sets (VCF & BED) for HG002, HG003, and HG004:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
```
Scripts assume these files are in `${TOP_DIR}/ref/GRCh37/`

## Software versions
* [`bcftools` 1.15.1](https://github.com/samtools/bcftools/releases/tag/1.15.1)

* [`bedtools` 2.30.0](https://github.com/arq5x/bedtools2/releases/tag/v2.30.0)

* [`samtools` 1.15.1](https://github.com/samtools/samtools/releases/tag/1.15.1)
