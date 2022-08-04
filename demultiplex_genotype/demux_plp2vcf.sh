#!/usr/bin/env bash
## Siva Kasinathan
# demux_plp2vcf.sh: Naively calls variants using samtools mpileup in the regions provided
# Usage: ./demux_plp2vcf.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     BAM: BAM-file to call variants
#     FASTA: GRCh37, hs37d5, reference genome fasta
#     CALLING_REGIONS: BED file containing shared set of GRCh37 verified calling regions produced by demux_process_truth_sets.sh
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/demultiplex_genotype/plp/:
#         BAM.q15.d2.vcf.gz: Naive VCF containing SNVs called in the provided regions with at least 2 reads and a minimum MAPQ of 15

TOPDIR=$1
BAM=$2
FASTA=$3
CALLING_REGIONS=$4

set -Eeuo pipefail

MINQ=15
MIND=2

PFX=`basename ${BAM}`

PLP_DIR=$TOP_DIR/analyses/HG/demultiplex_genotype/plp
mkdir -p $PLP_DIR

samtools mpileup \
    --no-BAQ \
    --fasta-ref $FASTA \
    --positions $CALLING_REGIONS \
    $BAM \
| bgzip \
> $PLP_DIR/${PFX%%.*}.plp.gz

zcat $PLP_DIR/${PFX%%.*}.plp.gz \
| $TOPDIR/../demultiplex_genotype/plp2vcf.py -q $MINQ -d $MIND - \
| bgzip \
> $PLP_DIR/${PFX%%.*}.q${MINQ}.d${MIND}.vcf.gz

bcftools index -t $PLP_DIR/${PFX%%.*}.q${MINQ}.d${MIND}.vcf.gz

echo "JOB DONE."