module load plink

HEADDIR=/oak/stanford/groups/smontgom/amarder/VariantPrioritization
mkdir -p $HEADDIR/1kg_maf
mkdir -p $HEADDIR/1kg_maf/eur
mkdir -p $HEADDIR/1kg_maf/all

chrNum=22

for chrNum in {1..22};
do

geno=/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels/chr$chrNum.eur.filtered
out=$HEADDIR/1kg_maf/eur/1kg_maf.chr$chrNum
 
plink --bfile $geno --freq --out $out

done


for chrNum in {1..22};
do

geno=/oak/stanford/groups/akundaje/soumyak/refs/plink_1kg_30x_hg38_snvs_indels/chr$chrNum.all.filtered
out=$HEADDIR/1kg_maf/all/1kg_maf.chr$chrNum

plink --bfile $geno --freq --out $out

done
