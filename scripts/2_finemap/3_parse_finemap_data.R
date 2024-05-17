trait <- "Lupus_Bentham_2015"

print("Loading packages and functions...")
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(susieR)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
source("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/finemap/finemap_pipeline_functions.R")
print("Packages and functions loaded.")

args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]

# stanford cluster-specific info
# can also hard call the GWAS file path and number of individuals in the GWAS
colocalRef <- fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv",data.table = F,stringsAsFactors = F)
colnames(colocalRef)[1] = "traitName" # compatibility issues
newGWAS <- as.logical(subset(colocalRef,traitName==trait)$newGWAS)
info.lst = marderstein_collect_info(trait,newGWAS)
filePath=info.lst[[1]]

################################################

finemap.save = list()

for (j in 1:2) {
# for (j in 1:1) { 
  if (j==1) {
    headdir = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_10/")
  } else if (j==2) {
    headdir = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_1/")
  } 
  flst = list.files(headdir)
  flst.cred = list.files()
  
  flst = flst[grep("variant_prob",flst)]
  f=sort(flst)
  
  f.cred = paste0(lapply(strsplit(f,"\\."),function(x) paste(x[1:2],collapse = ".")),".susie.cred_sets.txt")
  
  finemap.lst <- list()
  for (i in 1:length(f)) {
    print(paste0(i,"/",length(f)))
    f.iter = paste0(headdir,f[i])
    va = fread(f.iter,data.table = F,stringsAsFactors = F)
    f.iter = paste0(headdir,f.cred[i])
    cred = fread(f.iter,data.table = F,stringsAsFactors = F)
    if (j==2) {cred$Lval = 1} else if (!("Lval" %in% colnames(cred))) {cred$Lval = 10}
    print(cred[,c("sentinel","n_signals","n_snps_region","top_pip","converged","kb_window","Lval")])
    print(max(va$variable_prob))
    # if (is.na(cred$converged[1]) | cred$converged[1]) {
    finemap.lst[[i]] = va
    y = strsplit(f[i],"\\.")
    finemap.lst[[i]]$k <- lapply(y,function(x) x[[1]])
    finemap.lst[[i]]$sentinel <- lapply(y,function(x) x[[2]])
    finemap.lst[[i]]$converged <- cred$converged[1]
    finemap.lst[[i]]$Lval <- cred$Lval[1]
    # }
  }
  finemap = as.data.frame(do.call(rbind,finemap.lst))
  finemap[is.na(finemap$converged),"converged"] <- "NO_CS"
  finemap.save[[j]] = finemap
  
  # remove duplicates
  finemap.save[[j]] = (as.data.frame(as.data.table(finemap.save[[j]])[,.SD[which.max(variable_prob)],by='snp']))
}

finemap.save[[2]]$Lval=1
print("Merging L=1 and L=10...")
finemap.mg = merge(finemap.save[[1]],finemap.save[[2]],by='snp')
finemap.mg2 = finemap.mg[,c('sentinel.x','snp','converged.x','cs.x','variable_prob.x','converged.y','cs.y','variable_prob.y','Lval.x')]
colnames(finemap.mg2) = c('sentinel','snp','converged.L_10','CS.L_10','PIP.L_10','converged.L_1','CS.L_1','PIP.L_1','Lval')

####################

print("Reading LD expanded data...")
f.ldbud = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/all_ld_hg38/annot/",trait,".all_ld.r2.gene.txt")
df.ld_expand = fread(f.ldbud,data.table = F,stringsAsFactors = F)

####################

# path_to_downloaded_gwas="/oak/stanford/groups/smontgom/amarder/LDSC_pipeline/gwas/munge"
# filePath <- paste0(path_to_downloaded_gwas,'/hg38/',trait,'/',trait,'.txt.gz')
print("Reading in GWAS data...")
df.gwas <- fread(filePath,data.table = F,stringsAsFactors = F)
print("GWAS data loaded, remove duplicates from GWAS file...")
df.gwas = (as.data.frame(as.data.table(df.gwas)[,.SD[which.min(pvalue)],by='rsid']))

print("Merging GWAS and finemapping...")
df.mg = merge(df.gwas,finemap.mg2,by.x=c("rsid"),by.y=c("snp"),all.y=TRUE)

# sum(duplicated(finemap.mg2$snp))
# sum(duplicated(df.gwas$rsid))
# sum(duplicated(df.mg$rsid))

print("Reordering merged file by chr and SNP position...")
df.mg = df.mg[order(as.numeric(df.mg$chr),df.mg$snp_pos,decreasing = F),]

print("Saving:")
# f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/",trait,".susie_all.txt")
f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/all/",trait,".susie_all.txt")
fwrite(df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

print("Alright, we're done here. Vamanos!")

# colocalRef=/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv
# tail -n +4 $colocalRef | while read -r line;
# do
# trait=`echo $line | cut -f1 -d,`
# echo $trait
# cp /oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/$trait/SUSIE/$trait.susie_all.txt /oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/all/$trait.susie_all.txt
# done

