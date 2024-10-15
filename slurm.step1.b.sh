#!/bin/sh

#Step1.b
DIR1=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/MEE/step1/output4.vcf
OUT1=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/MEE/step1/SIOP_FAME_annotated.101524.vcf

sbatch vep.annot.step.1.b.shh ${DIR1} ${OUT1}
