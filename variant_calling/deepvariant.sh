#!/usr/bin/env bash
## Siva Kasinathan
# deepvariant.sh: Call variants from PacBio HiFi reads using deepvariant
# Usage: ./deepvariant.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     BAM: PacBio BAM file to be 
#     OUTPREFIX: Prefix to append to output files
#     FASTA: Genome FASTA file used for alignment
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/variant_calling/deepvariant/$OUTPREFIX
#          $OUTDIR/$OUTPREFIX.deepvariant.vcf.gz: Deepvariant VCF

TOPDIR=$1
BAM=$2
OUTPREFIX=$3
FASTA=$4


set -eu

THREADS=32

OUTDIR=${TOPDIR}/analyses/HG/variant_calling/deepvariant/${OUTPREFIX}
mkdir -p ${OUTDIR} ${OUTDIR}/log

run_deepvariant \
    --model_type PACBIO \
    --num_shards ${THREADS} \
    --verbosity 0 \
    --logging_dir ${OUTDIR}/log \
    --reads ${BAM} \
    --ref ${FASTA} \
    --output_vcf ${OUTDIR}/${OUTPREFIX}.deepvariant.vcf.gz

mv ${OUTDIR}/log/make_examples.log ${OUTDIR}/${OUTPREFIX}.make_examples.log
mv ${OUTDIR}/log/call_variants.log ${OUTDIR}/${OUTPREFIX}.call_variants.log
mv ${OUTDIR}/log/postprocess_variants.log ${OUTDIR}/${OUTPREFIX}.postprocess_variants.log

rm -rf ${OUTDIR}/log