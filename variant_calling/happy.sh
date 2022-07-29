#!/usr/bin/env bash

# Benchmark vcf against truth set using hap.py

VCF=$1
TRUTH_VCF=$2
TRUTH_BED=$3

set -eu

THREADS=12
OUTDIR=`dirname $VCF`
OUTPREFIX=`basename ${VCF%.vcf.gz}`.happy
FASTA=${HOME}/smrt_tag/reference/GRCh37/hs37d5.fa
STRAT=${HOME}/smrt_tag/reference/GRCh37/stratification/v3.0-GRCh37-v4.2.1-stratifications.tsv

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