phylop=/oak/stanford/groups/smontgom/amarder/data/phylop/hg38.phyloP100way.bedGraph

for chrnum in {1..22};
do

echo $chrnum

grep "^chr${chrnum}[[:space:]]" $phylop > /oak/stanford/groups/smontgom/amarder/data/phylop/hg38.phyloP100way.chr${chrnum}.bedGraph

done


#############

# srun --account=smontgom --partition=batch --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 

# This section uses bedtools intersect to merge variants w/ phylop, and account for some variants missing in phylop file.
# there does seem to just be some variants missing in the UCSC phylop, which is verified by CADD too

variant=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chrALL.filter.score.bed

for chrnum in {1..22};
do

echo $chrnum

variant2=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chr$chrnum.filter.score.v2.bed

cat $variant | grep "^chr${chrnum}[[:space:]]" | cut -f1,2,3,6 > $variant2

phylop=/oak/stanford/groups/smontgom/amarder/data/phylop/hg38.phyloP100way.chr${chrnum}.bedGraph
outfile=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chr${chrnum}.filter.score.v2.phylop.bed
outfile2=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chr${chrnum}.filter.score.v2.missing.phylop.bed

# there does seem to just be some variants missing in the UCSC phylop, which is verified by CADD too
conda activate bedtools
bedtools intersect -a $variant2 -b $phylop -wa -wb > $outfile

bedtools intersect -a $variant2 -b $phylop -v > $outfile2
awk '{print $0 "\tNA"}' $outfile2 > tmp
rm $outfile2
mv tmp $outfile2

cat $outfile | cut -f1,2,3,4,8 > tmp
cat $outfile2 | cut -f1,2,3,4,5 >> tmp
rm $outfile
rm $outfile2
mv tmp $outfile

bedtools sort -i $outfile > tmp
mv tmp $outfile

conda deactivate

rm $variant2

done



###############

outfile2=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chrALL.filter.score.v2.phylop.bed
rm $outfile2

for chrnum in {1..22};
do

echo $chrnum

outfile=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/lt_0.001/chr${chrnum}.filter.score.v2.phylop.bed

# cat $outfile | cut -f1,3,4,8 >> $outfile2
# cat $outfile | cut -f1,3,4,5 >> $outfile2
cat $outfile >> $outfile2

done




