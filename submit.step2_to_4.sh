#!/bin/sh
#SBATCH --job-name=step2    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1               # Job memory request
#SBATCH --output=Step2_to_4_%j.log   # Standard output and error log
#SBATCH --partition=medium


DIR1=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step1

#python step2.py ${DIR1}/SIOP_FAME_annotated.101524.vcf
python step2.py ${DIR1}/FAME.TopSignif_ENSEMBL_withNEAREST.100724.vcf
python step3.py
python step4.py