#!/bin/sh
#SBATCH --job-name=annotation    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=50G --nodes=1               # Job memory request
#SBATCH --output=annotation.vcf_%j.log   # Standard output and error log
#SBATCH --partition=short,long,medium

#DIR1=$1
#OUT1=$2

export outputFileName=SIOP_FAME_annotated.101524.vcf
DIR1=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step1/output4.vcf
OUT1=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step1/${outputFileName}
echo $OUT1
echo $DIR1

vep --cache --dir_cache /data/genomes/VEP/VEP_110dot1/GRCh38/cache/ENSEMBL/ --input_file ${DIR1} --output_file ${OUT1} --force_overwrite --everything --offline --nearest symbol --mane --dir_plugins /data/genomes/VEP/VEP_110dot1/GRCh38/plugins/  --fasta /data/Segre_Lab/users/yluo/Liftover/data/GRCh38.p13.genome.fa
