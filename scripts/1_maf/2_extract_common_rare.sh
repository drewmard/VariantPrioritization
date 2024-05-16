HEADDIR=/oak/stanford/groups/smontgom/amarder/VariantPrioritization
mkdir -p $HEADDIR/out/1kg_variants
mkdir -p $HEADDIR/out/1kg_variants/gt_0.1
mkdir -p $HEADDIR/out/1kg_variants/lt_0.01
mkdir -p $HEADDIR/out/1kg_variants/gt_0.05
mkdir -p $HEADDIR/out/1kg_variants/lt_0.005

for chrNum in {22..1};
do

echo $chrNum

awk 'NR==1 || ($5 < 0.01 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/all/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/lt_0.01/chr$chrNum.all.txt
awk 'NR==1 || ($5 < 0.01 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/eur/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/lt_0.01/chr$chrNum.eur.txt

awk 'NR==1 || ($5 > 0.1 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/all/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/gt_0.1/chr$chrNum.all.txt
awk 'NR==1 || ($5 > 0.1 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/eur/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/gt_0.1/chr$chrNum.eur.txt

awk 'NR==1 || ($5 < 0.005 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/all/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/lt_0.005/chr$chrNum.all.txt
awk 'NR==1 || ($5 < 0.005 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/eur/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/lt_0.005/chr$chrNum.eur.txt

awk 'NR==1 || ($5 > 0.05 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/all/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/gt_0.05/chr$chrNum.all.txt
awk 'NR==1 || ($5 > 0.05 && length($3)==1 && length($4)==1)' /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_maf/eur/1kg_maf.chr$chrNum.frq | awk -v OFS='\t' '{print $2, $5}' > $HEADDIR/out/1kg_variants/gt_0.05/chr$chrNum.eur.txt

done