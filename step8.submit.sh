#!/bin/bash
#SBATCH --job-name=step8   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=Step8_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

export RegeniefileName='SIOP_FAME.ALL.101524_responder'

python step8.py --regenie './step7/'${RegeniefileName}'.regenie' --outFile './step8_to_10/'${RegeniefileName}'.tsv'

