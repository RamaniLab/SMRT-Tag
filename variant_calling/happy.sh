#!/usr/bin/env bash
## Siva Kasinathan
# happy.sh: Benchmark vcf against truth set using hap.py
# Usage: ./happy.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     VCF: Deepvariant VCF file
#     TRUTH_VCF: Truth VCF file from GIAB
#     TRUTH_BED: Truth BED file from GIAB
#     FASTA: Genome FASTA file used for alignment
#     STRAT: Genome stratifications for GRCh37 from GIAB - v3.0-GRCh37-v4.2.1-stratifications.tsv
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/variant_calling/deepvariant/$OUTPREFIX:
#         $OUTDIR/$OUTPREFIX.deepvariant.happy: hap.py output file benchmarking variant calling

TOPDIR=$1
VCF=$2
TRUTH_VCF=$3
TRUTH_BED=$4
FASTA=$5
STRAT=$6

set -eu

THREADS=12
OUTDIR=`dirname $VCF`
OUTPREFIX=`basename ${VCF%.vcf.gz}`.happy

export HGREF=${FASTA}

hap.py \
    -r ${FASTA} \
    -o ${OUTDIR}/${OUTPREFIX} \
    -f ${TRUTH_BED} \
    --threads ${THREADS} \
    --pass-only \
    --engine=vcfeval \
    --verbose \
    --logfile ${OUTDIR}/${OUTPREFIX}.log \
    --stratification ${STRAT} \
    ${TRUTH_VCF} ${VCF}