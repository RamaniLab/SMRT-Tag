#!/bin/bash

set -Eeuo pipefail

HGS=(HG002 HG003 HG004)
MINQ=15
MIND=2

TOPDIR=${HOME}/smrt_tag/analyses/smrt_tag/demux
STRAT_BED=${HOME}/smrt_tag/reference/GRCh37/stratification/union/GRCh37_notinalldifficultregions.bed.gz
INDIR=${TOPDIR}/plp
OUTDIR=${TOPDIR}/plp_cmp
BEDDIR=${OUTDIR}/depth${MIND}_bed
mkdir -p $OUTDIR $BEDDIR

for AH in ${HGS[@]}; do
    VCF=${INDIR}/${AH}.q${MINQ}.d${MIND}.vcf.gz

    echo $VCF

    # Create and index depth BED; filter for stratification regions of interest
    PERBASE_BED=${TOPDIR}/${AH}/mosdepth/${AH}.per-base.bed.gz
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
        UNIQUE_SITES=${TOPDIR}/benchmark/unique/unique.${CH}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
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