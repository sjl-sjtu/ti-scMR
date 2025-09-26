#!/bin/bash

#SBATCH --job-name=split
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# convert illumina report to ped/map file: PLINK INput Report Plug-in in GenomeStudio

awk -v OFS="\t" '{gsub(",.*","",$2)}1' GSA-24v2-0_A1_b150_rsids.txt > names_rsid.txt

/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --pedmap onek1k_genotype \
    --allow-extra-chr \
    --make-pgen \
    --sort-vars \
    --out onek1k_genotype
/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --pfile onek1k_genotype \
    --rm-dup force-first \
    --snps-only just-acgt \
    --geno 0.03 --maf 0.01 --hwe 0.000001 \
    --update-name names_rsid.txt \
    --recode vcf \
    --out onek1k_genotype_converted

bgzip onek1k_genotype_converted.vcf
bcftools sort onek1k_genotype_converted.vcf.gz -Oz -o onek1k_genotype_converted.vcf.gz
tabix -p vcf onek1k_genotype_converted.vcf.gz
for i in {1..22}
do
  bcftools view onek1k_genotype_converted.vcf.gz -r ${i} -o splited/genotype_chr${i}.vcf.gz -Oz
done
bcftools view onek1k_genotype_converted.vcf.gz -r X -o splited/genotype_chrX.vcf.gz -Oz

# check strand, imputation, download
# check strand: conform-gt
for i in {1..22}
do
  java -jar conform-gt.jar ref=ref_vcf/chr${i}.1kg.phase3.v5a.vcf.gz gt=splited/genotype_chr${i}.vcf.gz chrom=i out=harmosed/geneotype_chr${i}
done
# imputation: Michigan Imputation Server (on line service)

cd ./imputed
bcftools concat *.vcf.gz -o merged.vcf.gz
bcftools query -l merged.vcf.gz > id_list.txt
bcftools query -f '%ID\n' merged.vcf.gz > snp_list.txt

/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf merged.vcf.gz \
    --update-ids ID_convert.txt \
    --keep id_slt.txt \
    --export vcf \
    --out merged_slt
/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf merged_slt.vcf \
    --maf 0.05 \
    --geno 0.05 \
    --hwe 0.000001 \
    --export A \
    --out merged_slt
cut merged_slt.raw -f 1-2 > slt_id_order.txt
head -n 1 merged_slt.raw | tr '\t' '\n' > common_snp_slt.txt

bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' merged_slt.vcf > imputed_snp_info.txt
/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf merged_slt.vcf \
  --pca

awk '
{
    for (i=1; i<=NF; i++) {
        if (NR == 1) {
            arr[i] = $i;
        } else {
            arr[i] = arr[i] " " $i;
        }
    }
}
END {
    for (i=1; i<=NF; i++) {
        print arr[i];
    }
}' merged_slt.raw > merged_slt_transpose.txt
