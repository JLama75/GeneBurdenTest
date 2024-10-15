import os
import sys
import argparse
import numpy as np
import pandas as pd

plink3='/modules/ogi-mbc/software/plink3/alpha3/bin/plink3'
plink='/modules/EasyBuild/software/PLINK/1.9b-6.10/plink'

#Set up argument parser
parser = argparse.ArgumentParser(description="Prepare files for Gene-burden test")


# Add arguments
parser.add_argument('--vcf', required=True, help="Path to the post-VariantQC VCF file")
parser.add_argument('--keep', required=True, help="list of samples to keep")
parser.add_argument('--monoThreshold', required=True, help="threshold for filtering out monomorphic variants")
parser.add_argument('--phenoFile', required=True, help="path to phenotype file with cases and controls")

# Parse the arguments
args = parser.parse_args()

output_path='./step1/'
pheno=args.phenoFile

out0=os.path.join(output_path, "output0")
out1=os.path.join(output_path, "output1")
out2=os.path.join(output_path, "output2")
out3=os.path.join(output_path, "output3")
out4=os.path.join(output_path, "output4")

#!/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --vcf {args.vcf} --keep {args.keep} --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 2000 --double-id --make-bed --out out1")
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --vcf {args.vcf} --keep {args.keep} --not-chr MT --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 2000 --double-id --make-bed --out {out0}")
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out0} --extract bed0 /data/Segre_Lab/users/yluo/resource/hg38_v0_exome_calling_regions.v1.interval_list.bed --make-bed --out {out1}")
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out1} --missing --out {out1}")
#monoThreshold=0.00093985
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out1} --geno 0.02 --make-bed --out {out2}")
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out2} --maf {args.monoThreshold} --make-bed --out {out3}")

#output4 are the variant list showing maf < 1% in our dataset.
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out3} --max-maf 0.01 --double-id --make-bed --out {out4}")
os.system(f"/modules/ogi-mbc/software/plink3/alpha3/bin/plink3 --bfile {out4} --export vcf-4.2  id-paste=iid --out {out4}")

##Then run sc_annnotation
file_vcf=os.path.join(output_path,'output4.vcf')
file_vcf_out=os.path.join(output_path,'output4.tsv')
file_vcf_out1=os.path.join(output_path,'output4_GT.tsv')
file_vcf_out1_case=os.path.join(output_path,'output4_GT_case.tsv')
file_vcf_out1_control=os.path.join(output_path,'output4_GT_control.tsv')

with open(file_vcf, encoding="ISO-8859-1") as vcf_in, open(file_vcf_out, 'w') as vcf_out:
    for line in vcf_in:
        if line[:2] != '##':
            tem_chr, tem_pos, tem_id, tem_ref, tem_alt, tem_qual, tem_filter, tem_info, tem_format, tem_GT = line.rstrip().split('\t', 9)
            new_line='\t'.join([tem_chr, tem_pos, tem_id, tem_ref, tem_alt, tem_info, tem_format, tem_GT])
            print(new_line, file=vcf_out)
            
print("outputing the extracted GT from vcf to tsv file: \n")           
with open(file_vcf_out) as vcf_in, open(file_vcf_out1, 'w') as vcf_out:
    for line in vcf_in:
        tem_chr, tem_pos, tem_id, tem_ref, tem_alt, tem_info, tem_format, tem_GT = line.rstrip().split('\t', 7)
        new_line=tem_id+'\t'+tem_GT
        print(new_line, file=vcf_out)
        
df_phenotype=pd.read_csv(pheno, sep=' ')
df_phenotype = df_phenotype[['IID', 'Phenotype']]
df_phenotype.columns=['IID','Phenotype']
print("phenotype file... \n")
print(df_phenotype)

df_sample_bfile=pd.read_csv('./step1/output1.fam', sep='\t', header=None)
df_sample_bfile.columns=['FID','IID','Var1','Var2','Var3','Var4']
print("output1.fam file... \n")
print(df_sample_bfile)

df_phenotype_subset=pd.merge(df_phenotype, df_sample_bfile, on=['IID'], how='inner')#removing excluded from all analysis 
df_phenotype_subset=df_phenotype_subset.loc[:,['IID','Phenotype']]
print("phenotype data with subsetted samples ... \n")
print(df_phenotype_subset)

df_GT=pd.read_csv(file_vcf_out1, sep='\t')
print("extracted genotype data ... \n")
print(df_GT)

df_control=df_phenotype_subset[df_phenotype_subset['Phenotype']==0]
df_control=df_control.loc[:,['IID']]

df_case=df_phenotype_subset[df_phenotype_subset['Phenotype']==1]
df_case=df_case.loc[:,['IID']]

df_GT_cases=df_GT[df_GT.columns[df_GT.columns.isin(df_case['IID'])]]
df_GT_cases.insert(0, "ID", df_GT['ID'])
print("Genotype information on cases: \n")
print(df_GT_cases)

df_GT_controls=df_GT[df_GT.columns[df_GT.columns.isin(df_control['IID'])]]
df_GT_controls.insert(0, "ID", df_GT['ID'])
print("Genotype information on controls: \n")
print(df_GT_controls)

print("outputing the extracted GT for controls and cases: \n")
df_GT_controls.to_csv(file_vcf_out1_control, sep='\t', index=False)
df_GT_cases.to_csv(file_vcf_out1_case, sep='\t', index=False)

