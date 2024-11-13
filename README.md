Pipeline:

1. sbatch step1.a.submit.sh (make sure to input the vcf file, file with samples you want to keep from the vcf file, your phenotype file and the monomorphic filter threshold )

		→ python step1.a.py --vcf ${your_vcf} --keep ${IDs_to_keep} --phenoFile ${your_phenotype_file} --monoThreshold 0.00093985 
2. sbatch vep.annot.step.1.b.sh (User have to specify the name of the output annotated vcf file.) For example:

		→ export outputFileName=your_annotated.vcf

Requires to activate VEP conda environment before running scirpt

4. sbatch submit.2_to_4.sh (will run step2, step3 and step4, user has to input path to annotated vcf within the bash script) For example:
   
		→ python step2.py $path_to_annotatedVCF
		

Use jupyter notebook for step 5,7, and 7 to generate AAF, covariates and phenotype files for regenie .

5. sbatch regenie.7.b.sh

6. sbatch submit.step8_to_9.sh (User has to specify the name of the output file from Regenie ran in previous step) For example: 
   
    		 → export RegeniefileName='SIOP_FAME.ALL.101524.responder'
   
		 → python step8.py --regenie './step7/'${RegeniefileName}'.regenie' --outFile  './step7/'${RegeniefileName}'.tsv'
   
   This will automatically submit 4 slurm jobs for 4 masks in parallel #User should edit the path to annotated vcf file and regenie output file in the Step.9.sh script. For Example:
   
			 →export vcf='./step1/SIOP_FAME_annotated.101524.vcf'
   
                         →export RegeniefileName='SIOP_FAME.ALL.101524.responder'
   
Use jupyter notebook for the final step 10






