#!/bin/sh

pathToDownload=$2
cd $pathToDownload

chrNum=$1
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrNum}.recalibrated_variants.vcf.gz

