eur_ids="/oak/stanford/groups/smontgom/amarder/bin/1kg/1kg_eur/eur_ids" # limit to european individuals
# GENOMEDIR="/oak/stanford/groups/smontgom/amarder/bin/1kg/30x/plink"
# the_1kg_prefix = "1kg.all.hg38_30x."

##################################################################################################

marderstein_collect_info = function(phenotype,newGWAS=NULL) {
  # Marderstein/stanford cluster-specific function to 
  # 1) identify GWAS filepath, and
  # 2) number of individuals in the GWAS.
  # alternatively, can also hard call both the GWAS file path and the number of individuals in the GWAS.
  # newGWAS is optional parameter, since it can also be found by querying the colocalRef file.
  
  print("marderstein_collect_info...")
  
  colocalRef <- fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv",data.table = F,stringsAsFactors = F)
  if (is.null(newGWAS <- as.logical(subset(colocalRef,trait==phenotype)$newGWAS)))
    
    if (newGWAS %in% c("TRUE","FALSE")) {newGWAS <- as.logical(newGWAS)} else {stop("newGWAS is neither TRUE nor FALSE.")}
  TMPDIR=paste0("/oak/stanford/groups/smontgom/amarder/tmp/",phenotype)
  
  filePath=NULL
  if (!newGWAS) {
    print("Not a newly downloaded GWAS...")
    df <- fread("/oak/stanford/groups/smontgom/amarder/bin/gwas_selection/gwas_selection - AllMungedGWAS.csv",data.table = F,stringsAsFactors = F)
    # re-name if there are duplicated "trait" names
    i <- which(duplicated(df$trait) | duplicated(df$trait,fromLast = T))
    df$trait[i] <- paste0(df$trait[i],'_',df$study[i],'_',df$author[i],'_',df$year[i])
    
    df.sub <- subset(df,trait==phenotype)
    print(phenotype)
    
    filePath=df.sub$path[1]
    NINDIV=as.numeric(df.sub$n_participants)
  } else {
    print("Alert: newly downloaded GWAS!")
    path_to_downloaded_gwas="/oak/stanford/groups/smontgom/amarder/LDSC_pipeline/gwas/munge"
    filePath <- paste0(path_to_downloaded_gwas,'/hg38/',phenotype,'/',phenotype,'.txt.gz')
    print(phenotype)
    
    NINDIV <- as.numeric(subset(colocalRef,trait==phenotype)$N_participants)
    # system(paste0("echo ",filePath," > ",TMPDIR,"/path_to_sumstat"))
  }
  return(list(filePath,NINDIV))
}

create_relevant_directories <- function(trait) {
  suppressMessages(suppressWarnings(dir.create("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap")))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait))))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/SNPLIST"))))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/LD"))))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/SUSIE"))))
  # suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/ldbuddy.L_1"))))
  # suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/snps/finemap/",trait,"/SUSIE/ldbuddy.L_10"))))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/SUSIE/window.L_1"))))
  suppressMessages(suppressWarnings(dir.create(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/SUSIE/window.L_10"))))
}

# check whether this SNP should be run:
bad_region_function = function(chrNum,posNum,remove_regions=NULL) {
  in_MHC = chrNum == 6 & posNum  > 25500000 & posNum < 33500000
  in_remove_regions = FALSE
  return(in_MHC | in_remove_regions)
}

# read in 1kg data
read1000G = function(chrNum) {
  print("Reading 1000G genotype data...")
  # f.geno = paste0(GENOMEDIR,"/",the_1kg_prefix,"chr",chrNum)
  f.geno = paste0("/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels/chr",chrNum,".eur.filtered")
  df.geno = fread(paste0(f.geno,".bim"),data.table = F,stringsAsFactors = F)
  return(df.geno)
}

# subset gwas file using lead snp
subset_gwas_using_lead_snp = function(df.gwas,chrNum,posNum,kbthres_to_use=1500) {
  M = 999999999
  while (M > 30000) {
    df.gwas.sub = subset(df.gwas,
                         chr==chrNum &
                           snp_pos >= posNum - as.numeric(kbthres_to_use)*1000 &
                           snp_pos <= posNum + as.numeric(kbthres_to_use)*1000
    )
    M = nrow(df.gwas.sub)
    if (M > 30000) {kbthres_to_use = kbthres_to_use - 500}
  }
  return(list(df.gwas.sub,kbthres_to_use))
}

# subset gwas file using lead snp
subset_gwas_using_lead_snp_ld_buddy = function(df.gwas,df.ld_expand.sub) {
  df.gwas.sub = subset(df.gwas,
                       chr==substring(df.ld_expand.sub$chr[1],4) &
                         snp_pos %in% df.ld_expand.sub[,"1_based_pos"] 
  )
  return(list(df.gwas.sub,NA))
}


reorganize_gwas = function(df.gwas.sub,chrNum) {
  
  # Reading 1000G genotype data...
  df.geno = read1000G(chrNum)
  
  print(paste0("Calculating z-scores..."))
  if (!("beta" %in% colnames(df.gwas.sub))) {
    if ("or" %in% colnames(df.gwas.sub)) {
      df.gwas.sub$beta = log(df.gwas.sub$or)
    }
  } 
  if ("beta" %in% colnames(df.gwas.sub)) {
    df.gwas.sub$z = df.gwas.sub$beta/df.gwas.sub$se
  } else if ("effect" %in% colnames(df.gwas.sub)) {
    df.gwas.sub$z = df.gwas.sub$effect/df.gwas.sub$se
  }
  print(paste0("Removing NA z-scores..."))
  df.gwas.sub = subset(df.gwas.sub,!is.na(z))
  
  # ####################################################
  # # Next, account for strand flip + reverse alleles:
  # # When comparing sumstats (GWAS) to genotype ref (1000G),
  # # A/C in sumstats and C/A in the ref is fine,
  # # the sumstats just need to be reversed (multiplying Z-score by -1).
  # # However, A/C and A/G is not fine.
  # # A/G and T/C is the same, just needs to be flipped according to strand.
  # # Same with A/G and C/T, but need to flip AND reverse.
  # # This function handles all the above.
  # # However:
  # # We are going to keep SNPs on ambiguous strands (e.g. A/T, C/G),
  # # but maybe these should be removed?
  # 
  # flip_reverse <- function(ss,info_ref) {
  #   flip_strand <- function(allele) {
  #     dplyr::case_when(
  #       allele == "A" ~ "T",
  #       allele == "C" ~ "G",
  #       allele == "T" ~ "A",
  #       allele == "G" ~ "C",
  #       TRUE ~ NA_character_
  #     )
  #   }
  #   
  #   # ss2 is original sumstats
  #   ss2 <- ss
  #   
  #   # ss3 is reversed sumstats
  #   ss3 <- ss
  #   ss3$effect_allele <- ss$non_effect_allele
  #   ss3$non_effect_allele <- ss$effect_allele
  #   ss3$z <- -ss$z
  #   
  #   #ss3 becomes original + reversed sumstats
  #   ss3 <- rbind(ss2, ss3) #####
  #   rm(ss2)
  #   
  #   # for non-ambiguous snps, do strand flip:
  #   ss4 <- ss3[!((ss3$non_effect_allele == "A" & ss3$effect_allele == "T") |
  #                  (ss3$non_effect_allele == "T" & ss3$effect_allele == "A") |
  #                  (ss3$non_effect_allele == "G" & ss3$effect_allele == "C") |
  #                  (ss3$non_effect_allele == "C" & ss3$effect_allele == "G")),]
  #   ss4$effect_allele <- flip_strand(ss4$effect_allele)
  #   ss4$non_effect_allele <- flip_strand(ss4$non_effect_allele)
  #   
  #   #ss4 becomes original + reversed + strand flip original + strand flip reverse sumstats
  #   ss4 <- rbind(ss3, ss4) ######
  #   ss.matched <- merge(ss4, info_ref[,c("snpid_ref","POS","A1","A2")], 
  #                       by.x=c("snp_pos","effect_allele","non_effect_allele"),by.y=c("POS","A2","A1"), all = FALSE)
  #   return(ss.matched)
  # }
  # df.geno = df.geno[,c(2,4,5,6)]
  # colnames(df.geno) <- c("snpid_ref","POS","A1","A2")
  # df.gwas.sub2 = flip_reverse(df.gwas.sub,info_ref = df.geno)
  # 
  # # remove duplicated SNPs:
  # df.gwas.sub2 <- subset(df.gwas.sub2,!duplicated(snpid_ref))
  
  df.gwas.sub2 = subset(df.gwas.sub,snp_id %in% df.geno[,2])
  
  # order SNPs by basepair position:
  df.gwas.sub2 <- df.gwas.sub2[order(df.gwas.sub2$snp_pos),]
  return(df.gwas.sub2)
}

computeLD = function(trait,lead_snp,i,chrNum,ldcalc=TRUE,eurONLY=TRUE) {
  
  print("Saving SNPLIST...")
  
  SNPLIST = paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/SNPLIST/",i,".",lead_snp,".snp.txt")
  tmpout = as.data.frame(df.gwas.sub2$snp_id)
  fwrite(tmpout,SNPLIST,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
  
  print("Saving LDDATA...")
  # f.geno = paste0(GENOMEDIR,"/",the_1kg_prefix,"chr",chrNum)
  f.geno = paste0("/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels/chr",chrNum,".eur.filtered")
  LDDATA = paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/LD/",i,".",lead_snp,".ld.txt")
  if (eurONLY) {
    cmd = paste0("mamba run -n plink plink --bfile ",f.geno," --keep ",eur_ids," --r square --keep-allele-order --extract ",SNPLIST," --out ",LDDATA)
  } else {
    cmd = paste0("mamba run -n plink plink --bfile ",f.geno," --r square --keep-allele-order --extract ",SNPLIST," --out ",LDDATA)
  }
  # cmd = paste0("mamba run -n plink plink --bfile ",f.geno," --r square --keep-allele-order --extract ",SNPLIST," --out ",LDDATA)
  if (ldcalc) {system(cmd)}
  
  print("reading in ld data into R...")
  R = fread(paste0(LDDATA,".ld"),data.table = F,stringsAsFactors = F)
  
  return(R)
}

LD_modify = function(R,z_scores,lead_snp_ind) {
  
  # This functions handles NA values in the LD matrix.
  # Why do NA values pop up?
  # In the first case, you have SNPs with no variation (e.g. all 0s; rs2952025). 
  # Here, the LD across all SNPs is "NA".
  # In a second case, you have two SNPs with variation (e.g. not all 0s; rs565686354 & rs73678004).
  # But, at least one SNP has a missing genotype call (e.g. an NA value; rs565686354).
  # Variation at the second SNP only occurs at the missing genotypes at the first SNP.
  # e.g. all 0s at rs565686354 except for 2 individuals: one has a 1 and the other is missing.
  # (cont.) rs73678004 is all 0s except for 1 individual; 
  # (cont.) this individual is also the one with missing rs565686354 call;
  # (cont.) so you end up with 0/0 (501), 0/1 (1), and 1/NA (1) which results in LD r = NA
  # To handle this:
  # - we first flag all SNPs with no variation (case 1).
  # - then we flag all problem SNP-SNP combinations, 
  # - and randomly remove 1 SNP from the SNP-SNP set (case 2)
  # - finally, we aim to retain the lead snp in the analysis
  
  # initialize
  j = rep(TRUE,nrow(R))
  to_drop=c() 
  
  if(any(is.na(R))) {
    cat(paste("> Warning! LD matrix contains missing values. Please wait while these are removed...\n"))
    
    # First, we are going to flag case 1 variants
    tmp = is.na(R)
    tmp = rowSums(tmp)
    # tmp = apply(R,1,function(x) sum(is.na(x)))
    if (sum(tmp==nrow(R)) > 0) {
      j = tmp!=nrow(R)
      R[!j,] <- 0
      R[,!j] <- 0
      # R <- R[j,j]
      # z_scores = z_scores[j]
      cat(paste("  ", sum(!j), "case 1 variants to be removed from final analysis...\n"))
    }
    
    
    # Next we flag case 2 variants:
    
    # convert upper triangle to -1 for now, s.t. we only account for NA LD values at 1 SNP
    tmp = R
    # tmp[lead_snp_ind,1] = NA # test for the "dont remove lead snp section
    tmp[upper.tri(tmp, diag = FALSE)] = -1
    
    # get row numbers of rows containing missing value
    to_drop <- which(is.na(as.data.frame(tmp)), arr.ind=TRUE)[,'row'] %>% unique()
    
    # don't remove lead snp - remove the other snp instead
    if (lead_snp_ind %in% to_drop) {
      to_drop = to_drop[!(to_drop %in% lead_snp_ind)]
      to_drop = unique(c(to_drop,which(is.na(tmp[lead_snp_ind,]))))
      # to_drop = unique(c(to_drop,which(is.na(R[lead_snp_ind,]))))
    }
    cat(paste("  ", length(to_drop), "case 2 variants to be removed from final analysis...\n"))
    
    # drop rows from matrix and sumstats
    to_drop = sort(c(to_drop,which(!j)))
    if (length(to_drop) > 0) {
      R = R[-to_drop, -to_drop]
      z_scores = z_scores[-to_drop]
    }
    cat(paste("  ", length(to_drop), "total variants removed from final analysis\n"))
  } else {
    z_scores = df.gwas.sub2$z
  }
  
  # find l, which represents which rows to keep in the gwas dataframe
  l = seq(1,nrow(df.gwas.sub2))
  if (length(to_drop) > 0) {l = l[-to_drop]}
  
  cat(paste("  ", nrow(R), "variants kept for final analysis\n"))
  
  return(list(R,z_scores,l))
}

process_susie_results = function(df,rss) {
  
  cat(paste0("Parsing susieR results...\n"))
  
  va<-summary(rss)$vars
  va$snp<-df$rsid[va$variable]
  va$chr<-df$chr[va$variable]
  va$pos<-df$snp_pos[va$variable]
  va = va[,2:4]
  
  #################
  
  converge_status = rss$converged
  log_bayes = as.data.frame(rss$lbf_variable)
  pip =  as.data.frame(rss$alpha)
  colnames(log_bayes) = df$rsid
  colnames(pip) = df$rsid
  
  suppressMessages(suppressWarnings(library(tidyr)))
  suppressMessages(suppressWarnings(library(dplyr)))
  null_cred_set = tibble(cs = NA) %>% 
    mutate(variant_ids = NA) %>% 
    mutate(variant_log10bf = NA) %>%
    mutate(variant_pip = NA) %>% 
    mutate(cs_size = 0) %>% 
    mutate(cs_log10bf = NA) %>%  
    mutate(cs_avg_r2 = NA) %>% 
    mutate(cs_min_r2 = NA) %>% 
    mutate(top_log10bf = NA) %>% 
    mutate(top_pip = NA) %>% 
    mutate(converged = NA) %>% 
    mutate(n_snps_region = nrow(df)) %>% 
    select(n_snps_region, cs_index=cs, cs_size, cs_log10bf:cs_min_r2, top_log10bf, top_pip, 
           variant_ids:variant_pip, converged)
  null_cred_set = as.data.frame(null_cred_set)
  
  if( !is.null(summary(rss)$cs) ){
    cred_set_summary = as_tibble(summary(rss)$cs) %>% 
      rowwise() %>% 
      
      # use the variant indices ("variables" in susie output) 
      # to extract the variant IDs from the input data and collapse to a single field
      mutate(variant_ids = paste0(df[as.numeric(unlist(strsplit(variable,","))),]$rsid, collapse=", ")) %>% 
      
      # use the variant IDs in each row (each independent credible set)
      # to look up the log10bf and posterior prob of inclusion (PIP) for each SNP in the credible set;
      # round these values to 3 sigfigs, and collapse to a single field
      mutate(variant_log10bf = paste0(round(log_bayes[cs, unlist(strsplit(variant_ids,", "))],3), collapse=", ")) %>%
      mutate(variant_pip = paste0(signif(pip[cs, unlist(strsplit(variant_ids,", "))],3), collapse=", ")) %>% 
      
      # Extract the maximum SNP-level log10bf and PIP for each credible set
      mutate(top_log10bf = max(as.numeric(unlist(strsplit(variant_log10bf,", "))))) %>% 
      mutate(top_pip = max(as.numeric(unlist(strsplit(variant_pip,", "))))) %>% 
      
      # calculate the size of each credible set as the number of SNPs included
      mutate(cs_size = length(unlist(strsplit(variable,",")))) %>%       
      
      # indicate if the credible set is viable (i.e. if the regression converged)
      mutate(converged = converge_status) %>% 
      
      # get number of SNPs tested by susieR
      mutate(n_snps_region = nrow(df)) %>% 
      select(n_snps_region, cs_index=cs, cs_size, cs_log10bf:cs_min_r2, top_log10bf, 
             top_pip, variant_ids:variant_pip, converged) #%>% 
    
    # add warning
    # mutate(warning = "NA")
    
  } else {
    cat(paste("> Warning! SusieR generated no credible sets\n"))
    cred_set_summary = null_cred_set #%>% mutate(warning = "no susieR credible sets")
  }
  
  cred_sets = cred_set_summary %>%
    mutate(gwas = trait, sentinel = lead_snp) %>%
    mutate(n_signals = nrow(cred_set_summary)) %>%
    select(gwas:n_signals, n_snps_region:converged) %>%
    mutate(kb_window = kbthres_to_use)
  cred_sets = as.data.frame(cred_sets)
  
  tmp = subset(df.gwas,rsid==lead_snp)
  cred_sets$chr = tmp$chr
  cred_sets$pos = tmp$snp_pos
  cred_sets$chrpos = paste0(subset(tmp,rsid==lead_snp)$chr,"_",subset(tmp,rsid==lead_snp)$snp_pos)
  cred_sets$numsnptested = nrow(va)
  
  #################
  
  return(list(va,cred_sets))
}

removeLD <- function(phenotype,lead_snp,i) {
  trait = phenotype
  LDDATA = paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/",trait,"/LD/",i,".",lead_snp,".ld.txt")
  cmd = paste0("rm ",LDDATA,".ld")
  system(cmd)
  print("LD file removed.")
}






