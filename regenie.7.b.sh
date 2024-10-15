#!/bin/bash
#SBATCH --job-name=Regenie  # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=5G --nodes=1              # Job memory request
#SBATCH --output=regenie.fame.responder_101524_%j.log   # Standard output and error log
#SBATCH --partition=short,long,medium

module load regenie/2.0.1
module load plink2/alpha5.11

regenie \
    --step 2 \
    --ignore-pred \
    --bed /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step1/output4 \
    --bt \
    --phenoFile /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step5_to_6/FAME.Responder.pheno.csv \
    --covarFile /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step5_to_6/FAME.GWAS_cov_pca_100924.csv \
    --phenoCol responder \
    --covarColList Age,Sex,steroid.medication,hx.glaucoma,drops.during6m,PC{1:10},count_synonymous \
    --set-list /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step2_to_4/study_set.tsv \
    --anno-file /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step2_to_4/study_annot.tsv \
    --minMAC 0.5 \
    --mask-def study_Mask.tsv \
    --build-mask 'max' \
    --aaf-file /data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step5_to_6/AAF_list.tsv \
    --out ./step7/SIOP_FAME.ALL.101524
