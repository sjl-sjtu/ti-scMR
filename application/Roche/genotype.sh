

awk -v OFS="\t" '{gsub(",.*","",$2)}1' GSA-24v3-0_A1_b151_rsids.txt > names_rsid.txt

/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf RES0103_GSAv3+_anon.vcf \
    --update-sex sexinfo.txt \
    --rm-dup force-first \
    --snps-only just-acgt \
    --geno 0.01 --maf 0.01 --hwe 0.000001 \
    --update-name names_rsid.txt \
    --recode vcf \
    --out MS_genotype

bgzip MS_genotype.vcf
bcftools sort MS_genotype.vcf.gz -Oz -o MS_genotype.vcf.gz
tabix -p vcf MS_genotype.vcf.gz
for i in {1..22}
do
  bcftools view MS_genotype.vcf.gz -r ${i} -o splited_hg38/genotype_chr${i}.vcf.gz -Oz
done

cd splited_hg38
bcftools norm -d both genotype_chr19.vcf.gz -Oz -o genotype_chr19_fl.vcf.gz

for i in (1..22) 
do
  java -jar /lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/conform-gt.24May16.cee.jar \
    ref=/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/genotype_chr/ref_vcf/b37.vcf/chr$i.1kg.phase3.v5a.vcf.gz \
    gt=splited/genotype_chr$i.vcf.gz \
    chrom=$i \
    match=ID \
    out=fliped_hg37/chr$i
done

cd fliped_hg37

for i in {1..22}
do
  tabix -p vcf chr${i}.vcf.gz
  bcftools annotate --rename-chrs chr_name_conv.txt chr${i}.vcf.gz -Oz -o ../rename/chr${i}.vcf.gz
done


# download imputation
cd imp_hg38

module load bcftools/1.16-gcc-8.5.0-openblas
bcftools concat *.vcf.gz -o merged.vcf.gz
bcftools query -l merged.vcf.gz > id_list.txt
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' merged.vcf.gz > imputed_snp_info.txt

/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf merged.vcf.gz \
    --maf 0.05 \
    --geno 0.05 \
    --hwe 0.000001 \
    --export A \
    --out merged_slt

head -n 1 merged_slt.raw | tr '\t' '\n' > common_snp_slt.txt

/lustre/home/acct-clsyzs/clsyzs/SunJianle/plink2/plink2 --vcf merged.vcf.gz \
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
