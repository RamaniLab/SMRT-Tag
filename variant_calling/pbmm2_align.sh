#!/usr/bin/env bash
## Siva Kasinathan
# pbmm2_align.sh: Align PacBio unaligned BAM files to a reference genome and run mosdepth to calculate depth stats.
# Usage: ./pbmm2_align.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     BAM: PacBio BAM file to be aligned
#     OUTPREFIX: Prefix to append to output files
#     OUTDIR: Directory to write output
#     FASTA: Reference genome FASTA file
#     REGIONS: Region file for depth calculations - GRCh37_notinalllowmapandsegdupregions.bed.gz
#
## Outputs:
#     Outputs are written to $OUTDIR
#         $OUTDIR/$OUTPREFIX.aligned.bam: Aligned BAM file
#         $OUTDIR/mosdepth/$OUTPREFIX.mosdepth* : mosdepth files including (per-base, regions, thresholds.bed etc.)

TOPDIR=$1
BAM=$2
OUTPREFIX=$3
OUTDIR=$4
FASTA=$5
REGIONS=$6

set -eu

DTYPE=smrt_tag
THREADS=32

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
    ${OUTDIR}/${OUTPREFIX}.align.sorted.bam

mosdepth \
    --threads ${THREADS} \
    --use-median \
    --by ${REGIONS} \
    ${OUTDIR}/mosdepth/${OUTPREFIX} \
    ${OUTDIR}/${OUTPREFIX}.align.sorted.bam