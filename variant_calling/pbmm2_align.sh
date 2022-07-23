#!/bin/bash

# Align reads using pbmm2 and compute per-base depth
# Input is pacbio BAM file with unaligned reads

BAM=$1
OUTPREFIX=$2

set -eu

DTYPE=smrt_tag
THREADS=32
TOPDIR=${HOME}/smrt_tag/analyses/${DTYPE}
OUTDIR=${TOPDIR}/${OUTPREFIX}
FASTA=${HOME}/smrt_tag/reference/GRCh37/hs37d5.fa
REGIONS=${HOME}/smrt_tag/reference/GRCh37/stratification/GRCh37_notinalllowmapandsegdupregions.bed.gz

mkdir -p ${OUTDIR}

pbmm2 align \
    --log-level INFO \
    --log-file ${OUTDIR}/${OUTPREFIX}.pbmm2_align.log \
    --preset HiFi \
    --sort \
    --num-threads ${THREADS} \
    --sample ${OUTPREFIX} \
    ${FASTA} \
    ${BAM} \
    ${OUTDIR}/${OUTPREFIX}.aligned.bam

mosdepth \
    --threads ${THREADS} \
    --use-median \
    --by ${REGIONS} \
    ${OUTDIR}/mosdepth/${OUTPREFIX} \
    ${OUTDIR}/${OUTPREFIX}.aligned.bam