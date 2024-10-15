#!/bin/bash
#SBATCH --job-name=step9   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=step9.extractGT.count_101524_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

#User should update the path to vcf and study name of regenie summary statistic

list_annot_run=$1 #list_annot_run=list_annot_run2 
out_name=$2 #Run2_Count, Run3_Count .. so on

#export vcf='./step1/SIOP_FAME_annotated.101524.vcf'
export vcf='./step1/FAME.TopSignif_ENSEMBL_withNEAREST.100724.vcf'
export RegeniefileName='SIOP_FAME.ALL.101524_responder'
regenie_out='./step8_to_10/'${RegeniefileName}'.tsv'
output_path='./step8_to_10/'

echo $regenie_out
echo running for $list_annot_run
python step9.py --annotatedVCF ${vcf} --regenie_output ${regenie_out} --output_path ${output_path} --output_name ${out_name} --list_annot_run ${list_annot_run}