#!/bin/bash

set -Eeuo pipefail

TOPDIR=${HOME}/smrt_tag/analyses/smrt_tag/demux/benchmark

HG002_VCF=${TOPDIR}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
HG003_VCF=${TOPDIR}/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
HG004_VCF=${TOPDIR}/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz

HG002_BED=${TOPDIR}/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
HG003_BED=${TOPDIR}/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed
HG004_BED=${TOPDIR}/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed

# Make calling regions
mkdir -p ${TOPDIR}/calling_regions

cat $HG002_BED $HG003_BED $HG004_BED \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge \
| bgzip \
> ${TOPDIR}/calling_regions/HG002-4.calling_regions.bed.gz

tabix ${TOPDIR}/calling_regions/HG002-4.calling_regions.bed.gz

# Get regions that are unique to each individual
mkdir -p ${TOPDIR}/unique

bcftools isec \
    --threads 4 \
    -n~100 -w 1 \
    -c some \
    -Oz -o ${TOPDIR}/unique/unique.HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

bcftools isec \
    --threads 4 \
    -n~010 -w 2 \
    -c some \
    -Oz -o ${TOPDIR}/unique/unique.HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

bcftools isec \
    --threads 4 \
    -n~001 -w 3 \
    -c some \
    -Oz -o ${TOPDIR}/unique/unique.HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

for VCF in ${TOPDIR}/unique/*.vcf.gz; do
    bcftools index -t ${VCF}
done

echo "JOB DONE."