#!/bin/bash
#SBATCH --job-name=Step1   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=2G --nodes=1              # Job memory request
#SBATCH --output=Step1.101424_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

vcf='/data/Segre_Lab/users/jlama/WES_new.ALL_050824/variantQC/SIOP.variant_qc.093024.filtered.vcf.bgz'
keep='/data/Segre_Lab/users/jlama/WES_new.ALL_050824/Phenotype/FAME.Keep.removeExcludeFromAllanlaysis_sham_EAS.100724.tsv'
pheno='/data/Segre_Lab/users/jlama/WES_new.ALL_050824/Phenotype/Allcohorts.Finalphenotype.100124.tsv'

python step1.a.py --vcf ${vcf} --keep ${keep} --phenoFile ${pheno} --monoThreshold 0.00093985 

