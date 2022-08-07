#!/usr/bin/env bash
## Siva Kasinathan
# demux_process_truth_sets.sh: Process reference regions and variants from HG{002,003,004} to produce shared calling regions and private variants
# Usage: ./demux_process_truth_sets.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     HG002_VCF: HG002 GRCh37 benchmark variant calls
#     HG003_VCF: HG003 GRCh37 benchmark variant calls
#     HG004_VCF: HG004 GRCh37 benchmark variant calls
#     HG002_BED: HG002 GRCh37 verified variant calling regions
#     HG003_BED: HG003 GRCh37 verified variant calling regions
#     HG004_BED: HG004 GRCh37 verified variant calling regions
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/:
#         HG002-4.calling_regions.bed.gz: Shared set of GRCh37 verified calling regions
#
#     Outputs are written to ${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/unique/
#         unique.HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz: Unique variants for HG002
#         unique.HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz: Unique variants for HG003
#         unique.HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz: Unique variants for HG004

set -Eeuo pipefail

TOPDIR=$1
HG002_VCF=$2
HG003_VCF=$3
HG004_VCF=$4
HG002_BED=$5
HG003_BED=$6
HG004_BED=$7

# Make calling regions
CALLING_DIR=${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/
mkdir -p $CALLING_DIR

cat $HG002_BED $HG003_BED $HG004_BED \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools merge \
| bgzip \
> ${CALLING_DIR}/HG002-4.calling_regions.bed.gz

tabix ${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/HG002-4.calling_regions.bed.gz

# Get regions that are unique to each individual
UNIQ_DIR=${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/unique/
mkdir -p $UNIQ_DIR

bcftools isec \
    --threads 4 \
    -n~100 -w 1 \
    -c some \
    -Oz -o ${UNIQ_DIR}/unique.HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

bcftools isec \
    --threads 4 \
    -n~010 -w 2 \
    -c some \
    -Oz -o ${UNIQ_DIR}/unique.HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

bcftools isec \
    --threads 4 \
    -n~001 -w 3 \
    -c some \
    -Oz -o ${UNIQ_DIR}/unique.HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
    $HG002_VCF $HG003_VCF $HG004_VCF

for VCF in ${UNIQ_DIR}/*.vcf.gz; do
    bcftools index -t ${VCF}
done


echo "JOB DONE."