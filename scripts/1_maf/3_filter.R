# conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/HarmonizeGWAS

library(data.table)

chrNum=22

for (set in c("lt_0.005","gt_0.05","lt_0.01","gt_0.1")) {
  print(set)
# for (set in c("lt_0.01","gt_0.1")) {
  for (chrNum in 22:1) {
    print(chrNum)
    f=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set,"/chr",chrNum,".eur.txt")
    frq_eur = fread(f,data.table = F,stringsAsFactors = F)
    colnames(frq_eur)[2] = "MAF_EUR"
    
    f=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set,"/chr",chrNum,".all.txt")
    frq_all = fread(f,data.table = F,stringsAsFactors = F)
    colnames(frq_all)[2] = "MAF_ALL"
    
    frq = merge(frq_eur,frq_all,by="SNP")
    f.out=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set,"/chr",chrNum,".filter.txt")
    fwrite(frq,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
}


# set="lt_0.005"
# f=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set)
# 
# set="lt_0.01"
# chrNum=22
# f=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set,"/chr",chrNum,".filter.txt")
# frq = fread(f,data.table = F,stringsAsFactors = F)
# frq = subset(frq,MAF_EUR < 0.005 & MAF_ALL < 0.005)
# f.out=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/",set,"/chr",chrNum,".filter.txt")
# fwrite(frq,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
# 



