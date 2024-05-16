library(data.table)

phenotype = "Alzheimers_Bellenguez_2022"
# phenotype <- "Lupus_Bentham_2015"
# newGWAS <- "FALSE"
args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]

print(paste(phenotype))

TMPDIR=paste0("/oak/stanford/groups/smontgom/amarder/tmp/",phenotype)

filePath=NULL
path_to_downloaded_gwas="/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/munge"
filePath <- paste0(path_to_downloaded_gwas,'/hg38/',phenotype,'/',phenotype,'.v2.txt.gz')

print(phenotype)

print("Reading GWAS...")
gwas <- fread(filePath,data.table = F,stringsAsFactors = F)

gwas.sub = subset(gwas,pvalue < 5e-8 & chr %in% c(1:22))
window_to_use = 1000000

gwas.sub.leads = c()
for (chrNum in unique(gwas.sub$chr)) {
  print(paste("chrNum",chrNum))
  # for (chrNum in unique(gwas.sub$chr)[1]) {
  gwas.sub2 = subset(gwas.sub,chr==chrNum)
  for (i in 1:nrow(gwas.sub2)) {
    if (i==1) {
      lead_data = gwas.sub2[i,]
    } else if (i > 1) {
      if ((abs(gwas.sub2[i,"snp_pos"]) - lead_data[,"snp_pos"]) > window_to_use) {
        gwas.sub.leads = c(gwas.sub.leads,lead_data[,"snp_id"])
        lead_data = gwas.sub2[i,]
      } else if (gwas.sub2[i,"pvalue"] < lead_data[,"pvalue"]) {
        lead_data = gwas.sub2[i,]
      }
    } 
  }
  gwas.sub.leads = c(gwas.sub.leads,lead_data[,"snp_id"])
}

library(data.table)
window_to_use=1000000
dir_to_use = paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/leadsnp/window_",as.character(as.integer(window_to_use)))
suppressWarnings(dir.create(dir_to_use))
f.out = paste0(dir_to_use,"/",phenotype,".lead.txt")
fwrite(as.data.frame(gwas.sub.leads),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
# subset(gwas.sub2,rsid %in% gwas.sub.leads)





