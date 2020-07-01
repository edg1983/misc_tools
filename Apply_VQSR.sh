BUNDLE="/well/gel/HICF2/ref/GATK_bundle/hg38"
input=$1
output=$2

#Make sites only VCF (discard sample level info)
#/apps/well/gatk/4.1.4.0/gatk MakeSitesOnlyVcf -I HICF2_RareDisease_GATKBestPractises.vcf.gz -O cohort_sitesonly.vcf.gz

#Compute VQSLOD tranches for INDELs
gatk --java-options "-Xmx128g -Xms24g" VariantRecalibrator \
    -V cohort_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
    -mode INDEL \
    --max-gaussians 4 \
    -resource mills,known=false,training=true,truth=true,prior=12:$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource axiomPoly,known=false,training=true,truth=false,prior=10:$BUNDLE/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource dbsnp,known=true,training=false,truth=false,prior=2:$BUNDLE/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O cohort_indels.recal \
    --tranches-file cohort_indels.tranches

#Compute VQSLOD tranches for SNPs
gatk --java-options "-Xmx128g -Xms24g" VariantRecalibrator \
    -V cohort_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -resource hapmap,known=false,training=true,truth=true,prior=15:$BUNDLE/hapmap_3.3.hg38.vcf.gz \
    -resource omni,known=false,training=true,truth=true,prior=12:$BUNDLE/1000G_omni2.5.hg38.vcf.gz \
    -resource 1000G,known=false,training=true,truth=false,prior=10:$BUNDLE/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource dbsnp,known=true,training=false,truth=false,prior=7:$BUNDLE/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O cohort_snps.recal \
    --tranches-file cohort_snps.tranches

#Apply VQSR filter to INDELs
gatk --java-options "-Xmx64g -Xms16g" \
    ApplyVQSR \
    -V $input \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode INDEL \
    -O indel.recalibrated.vcf.gz

#Apply VQSR filter to SNPs
gatk --java-options "-Xmx64g -Xms16g" \
    ApplyVQSR \
    -V indel.recalibrated.vcf.gz \
    --recal-file cohort_snps.recal \
    --tranches-file cohort_snps.tranches \
    --truth-sensitivity-filter-level 99.5 \
    --create-output-variant-index true \
    -mode SNP \
    -O $output
