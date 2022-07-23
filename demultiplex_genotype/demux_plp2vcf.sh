#!/bin/bash

BAM=$1

set -Eeuo pipefail

MINQ=15
MIND=2

TOPDIR=${HOME}/smrt_tag/analyses/smrt_tag/demux
FASTA=${HOME}/smrt_tag/reference/GRCh37/hs37d5.fa
CALLING_REGIONS=${TOPDIR}/benchmark/calling_regions/HG002-4.calling_regions.bed.gz

PFX=`basename ${BAM}`

mkdir -p plp

samtools mpileup \
    --no-BAQ \
    --fasta-ref $FASTA \
    --positions $CALLING_REGIONS \
    $BAM \
| bgzip \
> plp/${PFX%%.*}.plp.gz

zcat plp/${PFX%%.*}.plp.gz \
| plp2vcf.py -q $MINQ -d $MIND - \
| bgzip \
> plp/${PFX%%.*}.q${MINQ}.d${MIND}.vcf.gz

bcftools index -t plp/${PFX%%.*}.q${MINQ}.d${MIND}.vcf.gz

echo "JOB DONE."