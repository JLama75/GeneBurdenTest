#!/bin/bash
#SBATCH --job-name=step10   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=Step10.2_101524_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

echo -e "running step10..."
#dir=$1
#trait=$2
dir='./step8_to_10'
trait='FAME.Responders'
python step10.py --Dir ${dir} --traitName ${trait} 
