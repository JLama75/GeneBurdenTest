#!/bin/bash
#SBATCH --job-name=step10   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=Step10_101524_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

echo -e "running step10..."
python step10.py
