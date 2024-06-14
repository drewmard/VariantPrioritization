# conda activate r

library(data.table)

# 738882 Alzheimers_Bellenguez_2022.susie_all.txt
f = "/oak/stanford/groups/smontgom/amarder/neuro-variants old/output/snps/finemap/all/Alzheimers_Bellenguez_2022.susie_all.txt"
old_fm = fread(f,data.table = F,stringsAsFactors = F)

# 457255 Alzheimers_Bellenguez_2022.susie_all.txt
f = "/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/all/Alzheimers_Bellenguez_2022.susie_all.txt"
new_fm = fread(f,data.table = F,stringsAsFactors = F)
y = strsplit(new_fm$snp,":")
new_fm$chr = substring(unlist(lapply(y,function(x) x[1])),4)
new_fm$snp_pos = unlist(lapply(y,function(x) x[2]))

df = merge(old_fm,new_fm,by=c("chr","snp_pos"),all=TRUE)
sum(duplicated(df$snp[!is.na(df$snp)]))

table(df$PIP.L_1.x > 0.2,df$PIP.L_1.y > 0.2)

sum(df$CS.L_1.x > 0,na.rm = T)
sum(df$CS.L_1.y > 0,na.rm = T)

sum(df$PIP.L_1.x > 0.2,na.rm = T)
sum(df$PIP.L_1.y > 0.2,na.rm = T)

sum(df$PIP.L_1.x > 0.5,na.rm = T)
sum(df$PIP.L_1.y > 0.5,na.rm = T)

sum(df$PIP.L_1.x > 0.9,na.rm = T)
sum(df$PIP.L_1.y > 0.9,na.rm = T)

df[df$PIP.L_1.x > 0.2 & is.na(df$PIP.L_1.y),] # variants lost in new fine map

sum(df$CS.L_10.x > 0,na.rm = T)
sum(df$CS.L_10.y > 0,na.rm = T)

sum(df$PIP.L_10.x > 0.2,na.rm = T)
sum(df$PIP.L_10.y > 0.2,na.rm = T)

sum(df$PIP.L_10.x > 0.5,na.rm = T)
sum(df$PIP.L_10.y > 0.5,na.rm = T)

sum(df$PIP.L_10.x > 0.9,na.rm = T)
sum(df$PIP.L_10.y > 0.9,na.rm = T)

df[df$PIP.L_10.x > 0.2 & is.na(df$PIP.L_10.y),] # variants lost in new fine map
df[which(df$PIP.L_10.y > 0.2 & (is.na(df$PIP.L_10.x) | df$PIP.L_10.x < 0.2)),] # variants gained in new fine map


