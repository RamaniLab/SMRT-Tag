#!/bin/bash

# Call variants from PacBio HiFi reads using deepvariant

BAM=$1
OUTPREFIX=$2

set -eu

THREADS=32
TOPDIR=${HOME}/smrt_tag/analyses/smrt_tag
FASTA=${HOME}/smrt_tag/reference/GRCh37/hs37d5.fa

OUTDIR=${TOPDIR}/${OUTPREFIX}
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