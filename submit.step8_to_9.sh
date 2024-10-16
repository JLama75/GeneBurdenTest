#!/bin/bash
#SBATCH --job-name=step8_to_9   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=Step8_to_9_101524_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

export RegeniefileName='SIOP_FAME.ALL.101524_responder' #Change to your Regenie file name

echo -e "running step8 ..."
python step8.py --regenie ${RegeniefileName} 

echo -e "running step9..."
#Enter --vcf path_to_your_annotated_vcf 
./slurm.step9.sh --vcf './step1/FAME.TopSignif_ENSEMBL_withNEAREST.100724.vcf' --RegeniefileName ${RegeniefileName} 
