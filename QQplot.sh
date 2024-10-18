#!/bin/bash
#SBATCH --job-name=QQplot   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlama@meei.harvard.edu     # Where to send mail    
#SBATCH --mem=1G --nodes=1              # Job memory request
#SBATCH --output=QQplot.101724_%j.log   # Standard output and error log
#SBATCH --partition=medium,long,short

echo running QQplot for MEE
DIR=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/MEE/step8_to_10
Rscript no-save QQplot.R --file1 ${DIR}"/Run1_merged_Count.tsv" --file2 ${DIR}"/Run2_merged_Count.tsv" --file3 ${DIR}"/Run3_merged_Count.tsv" --file4 ${DIR}"/Run4_merged_Count.tsv" --trait "MEE.responder" --output ${DIR} 

echo running QQplot for FAME-EUR
DIR=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_EUR/step8_to_10
Rscript no-save QQplot.R --file1 ${DIR}"/Run1_merged_Count.tsv" --file2 ${DIR}"/Run2_merged_Count.tsv" --file3 ${DIR}"/Run3_merged_Count.tsv" --file4 ${DIR}"/Run4_merged_Count.tsv" --trait "FAME_EUR.responder" --output ${DIR} 

echo running QQplot for MEE-EUR
DIR=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/MEE_EUR/step8_to_10
Rscript no-save QQplot.R --file1 ${DIR}"/Run1_merged_Count.tsv" --file2 ${DIR}"/Run2_merged_Count.tsv" --file3 ${DIR}"/Run3_merged_Count.tsv" --file4 ${DIR}"/Run4_merged_Count.tsv" --trait "MEE_EUR.responder" --output ${DIR} 

echo running QQplot for FAME
DIR=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME_updated/step8_to_10
Rscript no-save QQplot.R --file1 ${DIR}"/Run1_merged_Count.tsv" --file2 ${DIR}"/Run2_merged_Count.tsv" --file3 ${DIR}"/Run3_merged_Count.tsv" --file4 ${DIR}"/Run4_merged_Count.tsv" --trait "FAME.responder" --output ${DIR} 
