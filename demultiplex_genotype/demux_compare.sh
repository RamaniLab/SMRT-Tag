#!/usr/bin/env bash
## Siva Kasinathan
# demux_compare.sh: Intersect called variants with private variants to determine degree of shared variants 
# Usage: ./demux_compare.sh 
#
## Inputs:
#     TOPDIR: $TOP_DIR
#     BAM: BAM-file to call variants
#     REGIONS: Region file for calculations - GRCh37_notinalllowmapandsegdupregions.bed.gz
#
## Outputs:
#     Outputs are written to ${TOPDIR}/analyses/HG/demultiplex_genotype/plp_cmp/:
#         depth2_bed/HG{002,003,004}.d2.GRCh37_notinalldifficultregions.bed.gz: BED file indicating regions where at least 2 reads are present
#         HG{002,003,004}.q15.d2_vs_HG{002,003,004}_unique.vcf.gz: VCF file containing variants present in a sequenced sample and private benchmark variants.
#         HG{002,003,004}.q15.d2_vs_HG{002,003,004}_unique.stats: Stats for above VCF file
#         HG{002,003,004}.d2_vs_HG{002,003,004}.base.vcf.gz: VCF file containing shared variants in experimental regions (at least 2 reads are present)
#         HG{002,003,004}.d2_vs_HG{002,003,004}.base.stats: Stats for above VCF file

set -Eeuo pipefail

HGS=(HG002 HG003 HG004)
MINQ=15
MIND=2

TOPDIR=$1
STRAT_BED=$2

INDIR=$TOP_DIR/analyses/HG/demultiplex_genotype/plp
OUTDIR=$TOP_DIR/analyses/HG/demultiplex_genotype/plp_cmp
BEDDIR=${OUTDIR}/depth${MIND}_bed
mkdir -p $OUTDIR $BEDDIR

for AH in ${HGS[@]}; do
    VCF=${INDIR}/${AH}.q${MINQ}.d${MIND}.vcf.gz

    echo $VCF

    # Create and index depth BED; filter for stratification regions of interest
    PERBASE_BED=${TOPDIR}/analyses/HG/demultiplex_genotype/align/mosdepth/${AH}.per-base.bed.gz ## mosdepth output from aligment
    DEPTH_BED=${BEDDIR}/${AH}.d${MIND}.$(basename ${STRAT_BED%.bed*}).bed.gz

    zcat $PERBASE_BED \
    | awk -v D=${MIND} '{if ($4 >= D) print}' \
    | bedtools merge -i - \
    | bedtools intersect \
        -u -a - -b $STRAT_BED \
    | bgzip \
    > $DEPTH_BED

    tabix $DEPTH_BED

    # Loop through and intersect with unique sites; compute stats
    for CH in ${HGS[@]}; do
        UNIQUE_SITES=${TOPDIR}/analyses/HG/demultiplex_genotype/calling_regions/unique/unique.${CH}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
        ISEC_VCF=${OUTDIR}/${AH}.q${MINQ}.d${MIND}_vs_${CH}_unique.vcf.gz

        echo "> vs $CH"

        # Perform intersect
        bcftools isec \
            --threads 8 \
            -n =2 -w 1 \
            -c some \
            --regions-file $DEPTH_BED \
            -Oz -o $ISEC_VCF \
            $VCF \
            $UNIQUE_SITES 

        bcftools index \
            -t -f \
            --threads 8 \
            $ISEC_VCF
        
        bcftools stats \
            --threads 8 \
            $ISEC_VCF \
        > ${ISEC_VCF%.vcf.gz}.stats

        # Get 'denominator' for intersects
        BASE_VCF=${OUTDIR}/${AH}.d${MIND}_vs_${CH}.base.vcf.gz

        bcftools view \
            --threads 8 \
            --regions-file $DEPTH_BED \
            -Oz -o $BASE_VCF \
            --types snps \
            $UNIQUE_SITES
        
        bcftools index \
            -t -f \
            --threads 8 \
            $BASE_VCF
        
        bcftools stats \
            --threads 8 \
            $BASE_VCF \
        > ${BASE_VCF%.vcf.gz}.stats

    done

done

echo "JOB DONE."