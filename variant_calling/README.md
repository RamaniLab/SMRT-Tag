# HG002 small variant (SNV & indel) calling

## Scripts
Please see [`SMRT-Tag_processing.ipynb`](../notebooks/SMRT-Tag_processing.ipynb) for examples of how these scripts are used.

Scripts assume the following top-level directory `${TOP_DIR}/analyses/HG/variant_calling/` with external files at paths noted below.

* [`pbmm2_align.sh`](./pbmm2_align.sh): Alignment with pbmm2 and coverage calculation with mosdepth.

* [`downsample_bam.sh`](./downsample_bam.sh): Downsample alignments to specified depth. Requires calculation of fraction of BAM reads to use for downsampling, which was done by dividing the desired depth (3, 5, 10, 15X) by `mosdepth` computed median depth (`grep -w total_region *.mosdepth.summary.txt | cut -f4`). For example, for HG002 SMRT-Tag and downsampling to 10X, the fraction of reads is 0.895 (=10/11.17).

* [`deepvariant.sh`](./deepvariant.sh): Call small variants with deepvariant.

* [`happy.sh`](./happy.sh): Benchmark using `hap.py`. For example, for SMRT-Tag deepvariant VCF at `${TOP_DIR}/smrt_tag/analyses/HG/variant_calling/HG002/HG002.deepvariant.vcf.gz`.

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
Scripts assume uncompressed tarball is at `${TOP_DIR}/ref/GRCh37/`

GIAB v4.2.1 benchmark call set (VCF & BED) for HG002:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
```
The VCF and BED need to be passed explicitly as above, place at `${TOP_DIR}/ref/`

HG002 PacBio Sequel II ~11 kb WGS reads aligned to hs37d5 from GIAB:
```
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
```

## Software versions
* [`deepvariant` 1.4.0](https://github.com/google/deepvariant/releases/tag/v1.4.0)

* [`hap.py` 0.3.12](https://github.com/Illumina/hap.py/releases/tag/v0.3.12)

* [`mosdepth` 0.3.3](https://github.com/brentp/mosdepth/releases/tag/v0.3.3)

* [`pbmm2` 1.9.0](https://github.com/PacificBiosciences/pbmm2/releases/tag/v1.9.0)

* [`samtools` 1.15.1](https://github.com/samtools/samtools/releases/tag/1.15.1)