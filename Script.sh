#!/bin/bash

# ==============================================================================
# 1. Reads mapping
# ==============================================================================

## 1.1 LRS data
meryl count k=15 output merylDB ${hg38_ref_genome}
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

winnowmap -W repetitive_k15.txt -L --cs --MD -t ${threads} -ax map-ont \
  ${hg38_ref_genome} ${sample}.fq.gz > ${sample}.sam

## 1.2 SRS data
### 1.2.1 Zhongyuan cohort
bwa-mem2 mem -t ${threads} -B 4 -O 6 -E 1 -M -R ${read_group_tag} \
  ${hg38_ref_genome} ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz > ${sample}.sam

### 1.2.2 Independent Han Chinese cohort
bwa mem -t ${threads} -B 4 -O 6 -E 1 -M -R ${read_group_tag} \
  ${hg38_ref_genome} ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz > ${sample}.sam


# ==============================================================================
# 2. Dark regions detection
# Scripts (DarkRegionFinder.jar, combine_DRF_output.py, create_annotation_bed.sh,
# camo_gene_pipeline.ogs) are adapted from Ebbert et al., Genome Biology (2019):
# https://github.com/mebbert/Dark_and_Camouflaged_genes/tree/master/scripts
# ==============================================================================

## 2.1 Run DarkRegionFinder (DRF) per sample
min_depth=-1          # Ensure no bases are written to low_cov_bed
min_mapq_mass=-1      # Ensure every position is written to low_mapq_bed
mapq_thresh=9         # GATK recommended cutoff
low_mapq_bed="${drf_output_dir}/${sample}.min_depth_${min_depth}.min_mapq_mass_${min_mapq_mass}.mapq_thresh_${mapq_thresh}.dark.low_mapq.bed"

java -Xmx32G -jar DarkRegionFinder.jar \
  -i ${sample}.bam \
  --human-ref ${hg38_ref_genome} \
  --min-region-size 1 \
  --mapq-threshold ${mapq_thresh} \
  --min-mapq-mass ${min_mapq_mass} \
  --min-depth ${min_depth} \
  --region-exclusivity \
  --low-coverage-bed-output /dev/null \
  --low-mapq-bed-output ${low_mapq_bed} \
  --incomplete-bed-output /dev/null

## 2.2 Combine DRF outputs across samples
mkdir -p ${result_dir} ${tmp_dir}

combine() {
  local suffix=$1
  local low_depth_out=$2
  local low_mapq_out=$3
  local combine_input="${tmp_dir}/${suffix}.input"

  find ${tmp_dir} -name "${suffix}" > ${combine_input}
  touch ${low_depth_out} ${low_mapq_out}

  python combine_DRF_output.py \
    ${combine_input} \
    ${low_depth_out} \
    ${low_mapq_out}
}

first_drf=$(ls ${drf_output_dir}/*.dark.low_mapq.bed | head -1)
total_lines=$(wc -l < ${first_drf})
chunks=16
nline=$(( (total_lines + chunks - 1) / chunks ))

for drf in ${drf_output_dir}/*.dark.low_mapq.bed; do
  base=$(basename ${drf})
  sample=${base%%.*}
  sample_dir="${tmp_dir}/${sample}/"
  mkdir -p ${sample_dir}
  split -l ${nline} ${drf} ${sample_dir} &
done
wait

low_depth_dir="${tmp_dir}/low_depth"
low_mapq_dir="${tmp_dir}/low_mapq"
mkdir -p ${low_depth_dir} ${low_mapq_dir}

for split_file in ${sample_dir}/*; do
  suffix=$(basename ${split_file})
  combine ${suffix} "${low_depth_dir}/${suffix}" "${low_mapq_dir}/${suffix}" &
done
wait

low_depth_merged="${result_dir}/${cohort_prefix}.dark.low_depth.bed"
low_mapq_merged="${result_dir}/${cohort_prefix}.dark.low_mapq.bed"
header="chrom\tstart\tend\tavg_nMapQBelowThreshold\tavg_depth\tavg_percMapQBelowThreshold"

echo -e "${header}" > ${low_depth_merged}
echo -e "${header}" > ${low_mapq_merged}
cat ${low_depth_dir}/* >> ${low_depth_merged}
cat ${low_mapq_dir}/* >> ${low_mapq_merged}

rm -rf ${tmp_dir}

## 2.3 Prepare gene annotation BED
bash create_annotation_bed.sh ${ensembl_grch38_gff3} ${anno_output_dir}

## 2.4 Identify dark and camouflaged regions
bash camo_gene_pipeline.ogs \
  -d ${low_depth_merged} \
  -m ${low_mapq_merged} \
  -g ${hg38_ref_genome} \
  -a ${grch38_annotation_bed} \
  -s ${sequencing_platform} \
  -v hg38 \
  -t ${threads} \
  -r ${result_dir}


# ==============================================================================
# 3. SNV calling
# ==============================================================================

## 3.1 LRS data
run_pepper_margin_deepvariant call_variant \
  -b ${sample}.bam \
  -f ${hg38_ref_genome} \
  -o ${output_dir} \
  -p ${sample} \
  -t ${threads} \
  --ont_r9_guppy5_sup \
  --gvcf

glnexus_cli \
  --config DeepVariantWGS \
  --threads ${threads} \
  --list ${gvcf_file_list} > ${lrs_jointcalling}.bcf

whatshap genotype \
  --chromosome ${chr} \
  --reference ${hg38_ref_genome} \
  -o ${sample}.${chr}.regenotype.vcf \
  ${sample}.${chr}.jointcalling.vcf \
  ${sample}.bam

## 3.2 SRS data
java -jar picard.jar MarkDuplicates \
  --INPUT ${sample}.bam \
  --OUTPUT ${sample}.dedup.bam \
  --METRICS_FILE ${sample}.metric.txt \
  --ASSUME_SORT_ORDER coordinate \
  --VALIDATION_STRINGENCY SILENT \
  --CREATE_INDEX true

gatk ApplyBQSR \
  -R ${hg38_ref_genome} \
  -I ${sample}.dedup.bam \
  -bqsr ${recal_table} \
  -O ${sample}.recal.bam

gatk HaplotypeCaller \
  -ERC GVCF \
  --intervals ${chr} \
  -R ${hg38_ref_genome} \
  --dbsnp ${dbsnp_grch38_vcf} \
  -I ${sample}.recal.bam \
  -O ${sample}.${chr}.g.vcf.gz

gatk CombineGVCFs \
  -R ${hg38_ref_genome} \
  -V ${gvcf_file_list} \
  -O combined.${chr}.g.vcf.gz \
  --intervals ${chr}

gatk GenotypeGVCFs \
  -R ${hg38_ref_genome} \
  -V combined.${chr}.g.vcf.gz \
  -O genotyped.${chr}.vcf.gz \
  --dbsnp ${dbsnp_grch38_vcf} \
  --intervals ${chr}


# ==============================================================================
# 4. Positive selection analysis
# ==============================================================================

## 4.1 Read-backed phasing (WhatsHap)
whatshap phase \
  -o ${sample}.phased.vcf \
  --reference=${hg38_ref_genome} \
  ${sample_unphased_vcf} ${sample}.bam

bgzip ${sample}.phased.vcf
tabix -f -p vcf ${sample}.phased.vcf.gz

bcftools merge \
  --threads ${threads} \
  -o ${pop}.whatshap.phased.vcf \
  -l ${whatshap_vcf_list}

vcftools --vcf ${pop}.whatshap.phased.vcf --chr ${chr} --recode --out ${pop}.whatshap.${chr}
bgzip ${pop}.whatshap.${chr}.recode.vcf
tabix -p vcf ${pop}.whatshap.${chr}.recode.vcf.gz

## 4.2 Statistical phasing (SHAPEIT4)
vcftools --gzvcf ${pop_unphased_vcf_gz} --chr ${chr} --recode --out ${pop}.unphased.${chr}
bgzip ${pop}.unphased.${chr}.recode.vcf
tabix -p vcf ${pop}.unphased.${chr}.recode.vcf.gz

shapeit4 \
  --input ${pop}.unphased.${chr}.recode.vcf.gz \
  --map ${shapeit4_b38_genetic_map} \
  --region ${chr} \
  --sequencing \
  --scaffold ${pop}.whatshap.${chr}.recode.vcf.gz \
  --output ${pop}.${chr}.phased.vcf.gz

tabix -p vcf -f ${pop}.${chr}.phased.vcf.gz

## 4.3 Genetic map interpolation (PLINK)
plink --vcf ${pop}.${chr}.phased.vcf.gz --make-bed --out ${pop}.${chr}
plink --bfile ${pop}.${chr} --cm-map ${plink_grch38_genetic_map} ${chr_num} --make-bed --out ${pop}.${chr}.withcm
cut -f 1-4 ${pop}.${chr}.withcm.bim > ${pop}.${chr}.withcm.map

## 4.4 iHS calculation and normalization (selscan)
selscan --ihs \
  --threads ${threads} \
  --vcf ${pop}.${chr}.phased.vcf.gz \
  --map ${pop}.${chr}.withcm.map \
  --out ${pop}.${chr}

norm --ihs \
  --files ${pop}.${chr}.ihs.out \
  --bins 100 \
  --min-snps 20 \
  --winsize 10000 \
  --bp-win


# ==============================================================================
# 5. Balancing selection analysis
# ==============================================================================

## 5.1 Input file preparation (glactools)
glactools vcfm2acf \
  --onlyGT \
  --fai ${hg38_ref_genome_fai} \
  ${pop_vcf} > ${pop}.${chr}.acf.gz

glactools acf2betascan \
  --fold ${pop}.${chr}.acf.gz \
  | bgzip -c > ${pop}.${chr}.beta.txt.gz

## 5.2 Balancing selection scan (BetaScan)
python BetaScan.py \
  -i ${pop}.${chr}.beta.txt.gz \
  -fold \
  -o ${pop}.${chr}.betascores.txt


# ==============================================================================
# 6. Construction of LRS-based reference panel
# ==============================================================================

## 6.1 Reference panel phasing (WhatsHap + SHAPEIT4)
# Performed per autosome (chr1-chr22) following the same workflow as Section 4.1-4.2.

## 6.2 Convert phased reference panel to M3VCF format (Minimac3)
Minimac3 \
  --refHaps ${ref_panel}.${chr}.phased.vcf.gz \
  --processReference \
  --prefix ${ref_panel}.${chr}


# ==============================================================================
# 7. Imputation to the independent Han Chinese cohort
# ==============================================================================

## 7.1 Target cohort phasing (Eagle)
eagle \
  --bfile ${target_bfile}.${chr} \
  --geneticMapFile ${eagle_hg38_genetic_map} \
  --chrom ${chr_num} \
  --outPrefix ${target_phased}.${chr} \
  --numThreads ${threads}

## 7.2 Whole-autosome SNV imputation (Minimac4)
minimac4 \
  --refHaps ${ref_panel}.${chr}.m3vcf.gz \
  --haps ${target_phased}.${chr}.vcf.gz \
  --chr ${chr} \
  --format GT,DS \
  --allTypedSites \
  --prefix ${target_imputed}.${chr}

## 7.3 Post-imputation quality control (BCFtools)
bcftools view \
  -i 'R2>=0.75 && MAF>=0.01' \
  -O z \
  -o ${target_imputed_qc}.${chr}.vcf.gz \
  ${target_imputed}.${chr}.dose.vcf.gz


# ==============================================================================
# 8. PheWAS and conditional analysis
# ==============================================================================

## 8.1 PheWAS (PLINK2)
plink2 \
  --bfile ${target_imputed_qc_bfile} \
  --pheno ${phenotypes} \
  --covar ${covariates} \
  --glm \
  --quantile-normalize \
  --out ${phewas_out}

## 8.2 Conditional analysis (GCTA-COJO)
gcta \
  --bfile ${target_imputed_qc_bfile} \
  --chr ${chr_num} \
  --maf 0.01 \
  --cojo-file ${summary_statistics}.ma \
  --cojo-slct \
  --out ${cojo_out}.${chr}
