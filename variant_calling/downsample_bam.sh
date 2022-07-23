#!/bin/bash

# Downsample a BAM file by selecting a fraction of reads

BAM=$1
DEPTH=$2
FRAC=$3
DTYPE=smrt_tag

set -eu

SEED=0
PREFIX=`basename ${BAM%%.*}`
EXTRA_THREADS=7

TOPDIR=${HOME}/smrt_tag/analyses/${DTYPE}
OUTDIR=${TOPDIR}/${PREFIX}/${DEPTH}X/seed${SEED}
OUTPREFIX=${PREFIX}.${DEPTH}X.seed${SEED}
REGIONS=${HOME}/smrt_tag/reference/GRCh37/stratification/GRCh37_notinalllowmapandsegdupregions.bed.gz

mkdir -p ${OUTDIR}

samtools view \
    --threads ${EXTRA_THREADS} \
    --subsample ${FRAC} \
    --subsample-seed ${SEED} \
    --bam \
    --with-header \
    --write-index \
    --output "${OUTDIR}/${OUTPREFIX}.bam##idx##${OUTDIR}/${OUTPREFIX}.bam.bai" \
    ${BAM}

# Sanity check: run mosdepth again to to check that downsampled
# depth is as expected

mosdepth \
    --threads ${EXTRA_THREADS} \
    --use-median \
    --by ${REGIONS} \
    ${OUTDIR}/mosdepth/${OUTPREFIX} \
    ${OUTDIR}/${OUTPREFIX}.bam