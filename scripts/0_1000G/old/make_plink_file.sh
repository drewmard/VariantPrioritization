#!/bin/sh

module load plink

for chrNum in {10..22}; do

# print iterator
echo $chrNum

# initialize
vcfGeno=/oak/stanford/groups/smontgom/amarder/bin/1kg/30x/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrNum}.recalibrated_variants.vcf.gz
plinkGeno=/oak/stanford/groups/smontgom/amarder/bin/1kg/30x/plink_indel/1kg.all.hg38_30x.chr${chrNum}
eur_ids=/oak/stanford/groups/smontgom/amarder/bin/1kg/1kg_eur/eur_ids2

# make plink file
plink --vcf $vcfGeno --make-bed --keep-allele-order --out $plinkGeno

# convert vcf to plink bed format, preserving ref/alt identities present in the vcf
awk '{print $1"\t"$1"_"$4"_"$6"_"$5"\t"$3"\t"$4"\t"$5"\t"$6}' $plinkGeno.bim > /oak/stanford/groups/smontgom/amarder/tmp/tmp

# replace snp id "." with identifier based on chr pos ref alt
rm $plinkGeno.bim
mv /oak/stanford/groups/smontgom/amarder/tmp/tmp $plinkGeno.bim

done

