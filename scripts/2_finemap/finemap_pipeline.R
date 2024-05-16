# code originally adapted from Page Goddard and Mike Gloudesmans!

# srun --account=default --partition=interactive --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

# trait <- "Lupus_Bentham_2015"; newGWAS = TRUE

print("Loading packages and functions...")
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(susieR)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
source("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/scripts/2_finemap/finemap_pipeline_functions.R")
print("Packages and functions loaded.")

# parameters:
# trait <- "Lupus_Bentham_2015"; newGWAS = TRUE
# trait="Alzheimers_Bellenguez_2022"
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
startInd <- ifelse(length(args) < 2,1,as.numeric(args[2]))

configFileName = "/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/config/munge.config"
config = fromJSON(configFileName)
filePath = paste0(config$output_base_dir,"hg38","/",trait,"/",trait,".txt.gz")

df = fread(,data.table = F,stringsAsFactors = F)

colocalRef <- fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv",data.table = F,stringsAsFactors = F)
colnames(colocalRef)[1] = "traitName" # compatibility issues
newGWAS <- as.logical(subset(colocalRef,traitName==trait)$newGWAS)


##########################################################################

# stanford cluster-specific info
# can also hard call the GWAS file path and number of individuals in the GWAS
info.lst = marderstein_collect_info(trait,newGWAS)
filePath=info.lst[[1]]
NINDIV=info.lst[[2]]

##########################################################################

# initialize
create_relevant_directories(trait)

# Reading data:
print("Reading GWAS file...")
df.gwas <- fread(filePath,data.table = F,stringsAsFactors = F)

# print("Reading LD buddies file...")
# f.ld_expand = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/all_ld_hg38/annot/",trait,".all_ld.r2.gene.txt")
# df.ld_expand=fread(f.ld_expand,data.table = F,stringsAsFactors = F)
# lead_snp_lst = unique(df.ld_expand$lead) # if only want to do GWAS snps, do this: unique(subset(df.ld_expand,rsid==lead & info=="GWAS")$lead)

f.lead = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/leadsnp/window_1000000/",trait,".lead.txt")
lead_snp_lst = fread(f.lead,data.table = F,stringsAsFactors = F,header = F)[,1]
M = length(lead_snp_lst)
print(paste0(M," unique GWAS loci total..."))

##########################################################################

print(paste0("Finemapping, starting at index ",startInd,"..."))
startInd=1
for (i in startInd:M) {
  
  lead_snp = lead_snp_lst[i]
  
  # Extract data
  autosomal=(subset(df.gwas,rsid==lead_snp)$chr %in% seq(1,22))
  if (!autosomal) {
    print(paste0('Skipping GWAS region ',i,' with sentinel SNP ',lead_snp,"."))
    print(paste0("GWAS region is a non-autosomal SNP (or has a messed up chromosome name)."))
    next
  }
  
  ####
  f.cred_sets = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_10/",i,".",lead_snp,".susie.cred_sets.txt")
  if (file.exists(f.cred_sets)) {
    cred_sets = fread(f.cred_sets,data.table = F,stringsAsFactors = F)
    if (cred_sets$converged[1] != "FALSE" | is.na(cred_sets$converged[1])) {next}
  }
  
  # input: df.sub
  df.sub = subset(df.gwas,rsid==lead_snp)
  bad_region = bad_region_function(chrNum=df.sub[,"chr"],posNum=df.sub[,"snp_pos"],remove_regions=NULL)
  # df.sub = subset(df.ld_expand,rsid==lead_snp)
  # bad_region = bad_region_function(chrNum=substring(df.sub$chr,4),posNum=df.sub[,"1_based_pos"],remove_regions=NULL)
  if (bad_region) {
    print(paste0('Skipping GWAS region ',i,' with sentinel SNP ',lead_snp,"."))
    print(paste0("GWAS region is in an invalid fine-mapping region (e.g. MHC)."))
    next
  } 
  
  print(paste0('Running GWAS region ',i,' with sentinel SNP ',lead_snp,"."))
  # df.ld_expand.sub = subset(df.ld_expand,lead==lead_snp)
  
  ######################################################
  
  # Subset & modify GWAS:
  gwas_subset.lst = subset_gwas_using_lead_snp(df.gwas,chrNum=df.sub[,"chr"],posNum=df.sub[,"snp_pos"],kbthres_to_use=1000)
  # gwas_subset.lst = subset_gwas_using_lead_snp(df.gwas,chrNum=substring(df.sub$chr,4),posNum=df.sub[,"1_based_pos"],kbthres_to_use=1500)
  df.gwas.sub = gwas_subset.lst[[1]]; kbthres_to_use = gwas_subset.lst[[2]]
  # df.gwas.sub2 = reorganize_gwas(df.gwas.sub,chrNum=substring(df.sub$chr,4))
  df.gwas.sub2 = reorganize_gwas(df.gwas.sub,chrNum=df.sub$chr)

  # Compute LD:
  # R = computeLD(trait,lead_snp,i,chrNum=substring(df.sub$chr,4))
  R = computeLD(trait,lead_snp,i,chrNum=df.sub$chr)
  z_scores = df.gwas.sub2$z
  lead_snp_ind = which(df.gwas.sub2$rsid==lead_snp)
  LD_modify.list = LD_modify(R,z_scores,lead_snp_ind); R = LD_modify.list[[1]]; z_scores = LD_modify.list[[2]]; l = LD_modify.list[[3]]
  
  ######################################################
  
  # run susie,
  print("Running susie window.L_10")
  Lval=10; did_converge=FALSE
  while (!did_converge & !is.na(did_converge)) {
    print(paste0("Using: L = ",Lval))
    rss <- tryCatch(
      { 
        susie_rss(z_scores, as.matrix(R), n=NINDIV, L = Lval) 
      },
      error = function(e) {
        print("running with full reference instead...")
        # Compute LD, except use all individuals instead of just EUR:
        R = computeLD(trait,lead_snp,i,chrNum=df.sub$chr,eurONLY=FALSE)
        # R = computeLD(trait,lead_snp,i,chrNum=substring(df.sub$chr,4),eurONLY=FALSE)
        z_scores = df.gwas.sub2$z
        lead_snp_ind = which(df.gwas.sub2$rsid==lead_snp)
        LD_modify.list = LD_modify(R,z_scores,lead_snp_ind); R = LD_modify.list[[1]]; z_scores = LD_modify.list[[2]]; l = LD_modify.list[[3]]
        # run susie:
        rss <- susie_rss(z_scores, as.matrix(R), n=NINDIV, L = Lval)
        return(rss)
      }
    )
    # process susie results,
    processed_susie = process_susie_results(df=df.gwas.sub2[l,],rss)
    va =   processed_susie[[1]]
    cred_sets =   processed_susie[[2]]
    cred_sets$Lval = Lval
    did_converge = cred_sets$converged[1]
    Lval=Lval-1
  }
  
  cred_sets
  
  f.va = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_10/",i,".",lead_snp,".susie.variant_prob.txt")
  fwrite(va,f.va,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
  f.cred_sets = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_10/",i,".",lead_snp,".susie.cred_sets.txt")
  fwrite(cred_sets,f.cred_sets,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

  print("Running susie window.L_1")
  rss <- susie_rss(z_scores, as.matrix(R), n=NINDIV, L = 1)
  # process susie results,
  processed_susie = process_susie_results(df=df.gwas.sub2[l,],rss)
  va =   processed_susie[[1]]
  cred_sets =   processed_susie[[2]]
  cred_sets
  
  # and save susie results!
  f.va = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_1/",i,".",lead_snp,".susie.variant_prob.txt")
  fwrite(va,f.va,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
  f.cred_sets = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/window.L_1/",i,".",lead_snp,".susie.cred_sets.txt")
  fwrite(cred_sets,f.cred_sets,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
  removeLD(trait,lead_snp,i)
}

cat("\n")
print("~~~~~~~FINEMAPPING COMPLETE~~~~~~~")
cat("\n")

