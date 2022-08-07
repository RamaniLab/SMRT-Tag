#!/usr/bin/env bash
## Siva Kasinathan
# downsample_bam.sh: Downsample aligned BAM files produced by pbmm2_align.sh
# Usage: ./downsample_bam.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     BAM: PacBio BAM file to be aligned
#     DEPTH: Depth of downsampling
#     FRAC: Fraction provided to samtools for downsampling
#     REGIONS: Region file for depth calculations - GRCh37_notinalllowmapandsegdupregions.bed.gz
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/variant_calling/downsample/${PREFIX}/${DEPTH}X/seed${SEED}
#         $OUTDIR/${PREFIX}.${DEPTH}X.seed${SEED}.bam: Downsampled aligned BAM file
#         $OUTDIR/mosdepth/$OUTPREFIX.mosdepth* : mosdepth files including (per-base, regions, thresholds.bed etc.)

TOPDIR=$1
BAM=$2
DEPTH=$3
FRAC=$4
REGIONS=$5
DTYPE=smrt_tag

set -eu

SEED=0
PREFIX=`basename ${BAM%%.*}`
EXTRA_THREADS=7

OUTDIR=${TOPDIR}/analyses/HG/variant_calling/downsample/${PREFIX}/${DEPTH}X/seed${SEED}
OUTPREFIX=${PREFIX}.${DEPTH}X.seed${SEED}

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