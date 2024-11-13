Pipeline:

1. sbatch step1.a.submit.sh

		→ python step1.a.py --vcf ${vcf} --keep ${keep} --phenoFile ${pheno} --monoThreshold 0.00093985 
2. sbatch vep.annot.step.1.b.sh 
		#User have to specify the name of the output annotated vcf file. For example:
		export outputFileName=SIOP_FAME_annotated.101524.vcf
		#Requires to activate VEP conda environment before running scirpt

3. sbatch submit.2_to_4.sh
		→ python step2.py $path_to_annotatedVCF
		→ python step3.py
		→ python step4.py

Use jupyter notebook for step 5,7, and 7 to generate AAF, covariates and phenotype files for regenie .

4. sbatch regenie.7.b.sh

5. sbatch submit.step8_to_9.sh
    		 → export RegeniefileName='SIOP_FAME.ALL.101524.responder'
   
		 → python step8.py --regenie './step7/'${RegeniefileName}'.regenie' --outFile  './step7/'${RegeniefileName}'.tsv'
   
   sbatch submit (will submit 4 slurm jobs for 4 masks in parallel) #User should edit the path to annotated vcf file and regenie output file in the Step.9.sh script. For Example:
   
		 →export vcf='./step1/SIOP_FAME_annotated.101524.vcf'
   
                 →export RegeniefileName='SIOP_FAME.ALL.101524.responder'
   
Use jupyter notebook for the final step 10






