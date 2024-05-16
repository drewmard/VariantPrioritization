#!/usr/bash

# Initialize
HEADDIR=/oak/stanford/groups/smontgom/amarder/VariantPrioritization
set=gt_0.05

for set in "gt_0.05" "gt_0.1" "lt_0.01" "lt_0.005";
do

echo $set

# Create file:
chrNum=1
echo $chrNum
f=$HEADDIR/out/1kg_variants/$set/chr$chrNum.filter.txt
cat $f > $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt

# Cycle through other chr:
for chrNum in {2..22};
do
echo $chrNum
f=$HEADDIR/out/1kg_variants/$set/chr$chrNum.filter.txt
tail -n+2 $f >> $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt

done

# Creating a scoring file:
# cat $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt | awk -v OFS='\t' 'NR>1 {split($1, parts, ":"); chr=parts[1]; pos=parts[2]; a1=parts[3]; a2=parts[4]; print chr, pos-1, pos, a1, a2, $1}' $HEADDIR/out/1kg_variants/$set/chrALL.filter.score.bed
cat $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt | awk -v OFS='\t' 'NR>1 {split($1, parts, ":"); chr=parts[1]; pos=parts[2]; a1=parts[3]; a2=parts[4]; print chr, pos, a1, a2, $1}' > $HEADDIR/out/1kg_variants/$set/chrALL.filter.score_input.txt

echo "${HEADDIR}/out/1kg_variants/${set}/chrALL.filter.txt is done!"
done


mkdir -p $HEADDIR/out/1kg_variants/lt_0.001
awk 'NR==1 || ($2 < 0.001 && $3 < 0.001)' $HEADDIR/out/1kg_variants/lt_0.005/chrALL.filter.txt > $HEADDIR/out/1kg_variants/lt_0.001/chrALL.filter.txt
set=lt_0.001
# cat $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt | awk -v OFS='\t' 'NR>1 {split($1, parts, ":"); chr=parts[1]; pos=parts[2]; a1=parts[3]; a2=parts[4]; print chr, pos, a1, a2, $1}' > $HEADDIR/out/1kg_variants/$set/chrALL.filter.score_input.txt
chmod 777 $HEADDIR/out/1kg_variants/$set/chrALL.filter.score_input.txt
set=gt_0.05
# cat $HEADDIR/out/1kg_variants/$set/chrALL.filter.txt | awk -v OFS='\t' 'NR>1 {split($1, parts, ":"); chr=parts[1]; pos=parts[2]; a1=parts[3]; a2=parts[4]; print chr, pos, a1, a2, $1}' > $HEADDIR/out/1kg_variants/$set/chrALL.filter.score_input.txt
chmod 777 $HEADDIR/out/1kg_variants/$set/chrALL.filter.score_input.txt
