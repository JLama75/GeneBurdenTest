import os
import sys
import argparse
import numpy as np
import pandas as pd
import re
from collections import defaultdict
import json
import subprocess
import scipy
from collections import Counter
from io import StringIO
from scipy.stats import chi2_contingency

#Set up argument parser
parser = argparse.ArgumentParser(description="Extract GT for each Run")

#renaming output_skato to regenie_output
#regenie_output='../step9_to_11/SIOP_FAME.ALL.responder.100924_responder_Phenotype.tsv'
# Add arguments
parser.add_argument('--annotatedVCF', required=True, help="Path to the annotated VCF file")
parser.add_argument('--regenie_output', required=True, help="Path to the Regenie output file")
parser.add_argument('--output_path', required=True, help="Directory path where output files will be saved")
parser.add_argument('--output_name', required=True, help="Name for the output files")
parser.add_argument('--list_annot_run', required=True, help="Name of the annotation run out of four runs: list_annot_run1, list_annot_run2, list_annot_run3 or list_annot_run4")

#Add {args.annotatedVCF}
# Parse the arguments
args = parser.parse_args()

# Ensure output directory exists
if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

file_vcf=args.annotatedVCF

list_annot_run=args.list_annot_run

print(list_annot_run)
#list_annot_run1
print(args.output_name)
GT_control='./step1/output4_GT_control.tsv'
GT_case='./step1/output4_GT_case.tsv'

list_annot_directory='./step2_to_4/' #list_LOF1.tsv, list_LOF2.tsv ...so on
suffix='.tsv'

geneset_LoF1='./step2_to_4/study_LoF1.SetID'
geneset_LoF2='./step2_to_4/study_LoF2.SetID'
geneset_missense='./step2_to_4/study_missense.SetID'
geneset_moderate='./step2_to_4/study_moderate.SetID'
geneset_modifier='./step2_to_4/study_modifier.SetID'
geneset_low='./step2_to_4/study_low.SetID'
geneset_synonymous='./step2_to_4/study_synonymous.SetID'

list_annot_LoF1=['LoF1']
list_annot_LoF2=['LoF2']
list_annot_missense=['missense']
list_annot_moderate=['moderate']
list_annot_modifier=['modifier']
list_annot_low=['low']
list_annot_synonymous=['synonymous']

list_annot_run1=['LoF1']
list_annot_run2=['LoF1','LoF2','missense']
list_annot_run3=['LoF1','LoF2','missense','moderate','modifier','low']
list_annot_run4=['synonymous']

df_geneset_LoF1=pd.read_csv(geneset_LoF1, sep='\t', header=None)
df_geneset_LoF2=pd.read_csv(geneset_LoF2, sep='\t', header=None)
df_geneset_missense=pd.read_csv(geneset_missense, sep='\t', header=None)
df_geneset_moderate=pd.read_csv(geneset_moderate, sep='\t', header=None)
df_geneset_modifier=pd.read_csv(geneset_modifier, sep='\t', header=None)
df_geneset_low=pd.read_csv(geneset_low, sep='\t', header=None)
df_geneset_synonymous=pd.read_csv(geneset_synonymous, sep='\t', header=None)


df_geneset_LoF1.columns=['GENE','SNP']
df_geneset_LoF2.columns=['GENE','SNP']
df_geneset_missense.columns=['GENE','SNP']
df_geneset_moderate.columns=['GENE','SNP']
df_geneset_modifier.columns=['GENE','SNP']
df_geneset_low.columns=['GENE','SNP']
df_geneset_synonymous.columns=['GENE','SNP']
print("df_Geneset_LOF1: \n")
print(df_geneset_LoF1)

list_geneset_LoF1=df_geneset_LoF1['SNP'].tolist()
list_geneset_LoF2=df_geneset_LoF2['SNP'].tolist()
list_geneset_missense=df_geneset_missense['SNP'].tolist()
list_geneset_moderate=df_geneset_moderate['SNP'].tolist()
list_geneset_modifier=df_geneset_modifier['SNP'].tolist()
list_geneset_low=df_geneset_low['SNP'].tolist()
list_geneset_synonymous=df_geneset_synonymous['SNP'].tolist()

dicts_GT_count_geneset={}
dicts_GT_count_geneset['LoF1']=list_geneset_LoF1
dicts_GT_count_geneset['LoF2']=list_geneset_LoF2
dicts_GT_count_geneset['missense']=list_geneset_missense
dicts_GT_count_geneset['moderate']=list_geneset_moderate
dicts_GT_count_geneset['modifier']=list_geneset_modifier
dicts_GT_count_geneset['low']=list_geneset_low
dicts_GT_count_geneset['synonymous']=list_geneset_synonymous

#list_annot_run=list_annot_run2
#list_annot_run=[list_annot_run2, list_annot_run3, list_annot_run4]

#Define different categories
list_LoF1=['stop_gained','splice_acceptor_variant','splice_donor_variant','frameshift_variant','transcript_ablation']
list_LoF2=['stop_lost','start_lost','transcript_amplification','feature_elongation','feature_truncation']
list_missense=['missense_variant']
list_moderate=['inframe_insertion','inframe_deletion','protein_altering_variant']
list_modifier=['NMD_transcript_variant','3_prime_UTR_variant','5_prime_UTR_variant','TF_binding_site_variant','intergenic_variant','intron_variant','mature_miRNA_variant','non_coding_transcript_exon_variant','non_coding_transcript_variant','downstream_gene_variant','regulatory_region_variant','upstream_gene_variant','coding_sequence_variant','TFBS_ablation','TFBS_amplification','regulatory_region_ablation','regulatory_region_amplification','intergenic_variant','sequence_variant']
list_low=['splice_donor_5th_base_variant','splice_region_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant','incomplete_terminal_codon_variant']
list_synonymous=['synonymous_variant','start_retained_variant','stop_retained_variant']
dicts_annotation={'MANE':0,'LoF1':1,'LoF2':2,'missense':3,'moderate':4,'modifier':5,'low':6,'synonymous':7}

list_LoF1=['MANE_' + x for x in list_LoF1 if str(x)] + list_LoF1
list_LoF2=['MANE_' + x for x in list_LoF2 if str(x)] + list_LoF2
list_missense=['MANE_' + x for x in list_missense if str(x)] + list_missense
list_moderate=['MANE_' + x for x in list_moderate if str(x)] + list_moderate
list_modifier=['MANE_' + x for x in list_modifier if str(x)] + list_modifier
list_low=['MANE_' + x for x in list_low if str(x)] + list_low
list_synonymous=['MANE_' + x for x in list_synonymous if str(x)] + list_synonymous

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

summary_GT_control=pd.DataFrame()
with open(GT_control) as inny:
    list_control_GT=[]
    list_control_nonMiss=[]
    list_ID=[]
    list_individual=[]
    for line in inny:
        temp_id=line.rstrip().split('\t')[0]
        temp_GT=line.rstrip().split('\t')[1:]
        if temp_id=='ID':
            list_samples=temp_GT
        else:
            index=[]
            if '0/1' in temp_GT:
                index.append(temp_GT.index('0/1'))
                #display(index)
            if '1/0' in temp_GT:
                index.append(temp_GT.index('1/0'))
                #display(index)
            if '1/1' in temp_GT:
                index.append(temp_GT.index('1/1'))
                #display(index)
            #display(index)
            list_individual_curr = [list_samples[i] for i in index]
            list_individual.append(list_individual_curr)
            if temp_GT.count('0/0') < temp_GT.count('1/1'):
                count_GT=temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('0/0')*2
            elif temp_GT.count('0/0') >= temp_GT.count('1/1'): 
                count_GT=temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('1/1')*2
            count_nonMiss=2*(temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('1/1')+temp_GT.count('0/0'))
            list_control_GT.append(count_GT)
            list_control_nonMiss.append(count_nonMiss)
            list_ID.append(temp_id)
summary_GT_control['SNP']=list_ID
summary_GT_control['CountGT_control']=list_control_GT
summary_GT_control['nonMiss_control']=list_control_nonMiss
summary_GT_control['individual_control']=list_individual


summary_GT_case=pd.DataFrame()
with open(GT_case) as inny:
    list_case_GT=[]
    list_case_nonMiss=[]
    list_ID=[]
    list_individual=[]
    for line in inny:
        temp_id=line.rstrip().split('\t')[0]
        temp_GT=line.rstrip().split('\t')[1:]
        if temp_id=='ID':
            list_samples=temp_GT
        else:
            index=[]
            if '0/1' in temp_GT:
                index.append(temp_GT.index('0/1'))
            if '1/0' in temp_GT:
                index.append(temp_GT.index('1/0'))
            if '1/1' in temp_GT:
                index.append(temp_GT.index('1/1'))
            list_individual_curr = [list_samples[i] for i in index]
            list_individual.append(list_individual_curr)
            if temp_GT.count('0/0') < temp_GT.count('1/1'):
                count_GT=temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('0/0')*2
            elif temp_GT.count('0/0') >= temp_GT.count('1/1'):
                count_GT=temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('1/1')*2
            count_nonMiss=2*(temp_GT.count('0/1')+temp_GT.count('1/0')+temp_GT.count('1/1')+temp_GT.count('0/0'))
            list_case_GT.append(count_GT)
            list_case_nonMiss.append(count_nonMiss)
            list_ID.append(temp_id)
summary_GT_case['SNP']=list_ID
summary_GT_case['CountGT_case']=list_case_GT
summary_GT_case['nonMiss_case']=list_case_nonMiss
summary_GT_control['individual_case']=list_individual

summary_GT_case_control=pd.merge(summary_GT_control,summary_GT_case,on=['SNP'],how='inner')

print("summary_GT_case_control dataframe: \n")
print(summary_GT_case_control)

dicts=nested_dict(2,list)
with open(file_vcf) as vcf_in:
    for line in vcf_in:
        if line[0]=='#':
            continue
        tmp_Variation,tmp_Location,tmp_Allele,tmp_Gene,tmp_Feature,tmp_Feature_type,tmp_Consequence,tmp_cDNA_position,tmp_CDS_position,tmp_Protein_position,tmp_Amino_acids,tmp_Codons,tmp_Existing_variation,tmp_Extra = line.rstrip().split('\t')
        tmp_annot=tmp_Consequence
        if 'SYMBOL=' in tmp_Extra:
            result=re.search('SYMBOL=(.*);SYMBOL_SOURCE', tmp_Extra)
            tmp_gene=result.group(1)
        #if more than one annotation,should be assign the most deleterious one
        list_annotation=tmp_annot.split(',')
        if 'MANE' not in tmp_Extra:
            if len(list_annotation) == 1:
                dicts[tmp_Variation][tmp_gene].append(list_annotation[0])
            elif len(list_annotation) > 1:
                for annot_in_list in list_annotation:
                    dicts[tmp_Variation][tmp_gene].append(annot_in_list)
        elif 'MANE' in tmp_Extra:
            if len(list_annotation) == 1:
                dicts[tmp_Variation][tmp_gene].append('MANE_'+list_annotation[0])
            elif len(list_annotation) > 1:
                for annot_in_list in list_annotation:
                    dicts[tmp_Variation][tmp_gene].append('MANE_'+annot_in_list)

for key, val in dicts.items():
    for val_key, val_val in val.items():
        list_tmp=[]
        flag_mane=0
        for element_annotation in val_val:
            if 'MANE' in element_annotation and flag_mane==0:
                list_tmp=[]
                list_tmp.append(element_annotation)
                flag_mane=1
            elif 'MANE' in element_annotation and flag_mane==1:
                list_tmp.append(element_annotation)
                flag_mane=1
            elif 'MANE' not in element_annotation and flag_mane==0:
                list_tmp.append(element_annotation)
        dicts[key][val_key]=list_tmp

for key, val in dicts.items():
    for val_key, val_val in val.items():
        tmp_deleterious=''
        for element_annotation in val_val:
            if element_annotation in list_LoF1:
                element_annotation_category = 'LoF1'
            elif element_annotation in list_LoF2:
                element_annotation_category = 'LoF2'
            elif element_annotation in list_missense:
                element_annotation_category = 'missense'
            elif element_annotation in list_moderate:
                element_annotation_category = 'moderate'
            elif element_annotation in list_modifier:
                element_annotation_category = 'modifier'
            elif element_annotation in list_low:
                element_annotation_category = 'low'
            elif element_annotation in list_synonymous:
                element_annotation_category = 'synonymous'
            if tmp_deleterious == '':
                tmp_deleterious=element_annotation_category
            elif tmp_deleterious != '':
                if dicts_annotation[element_annotation_category] <= dicts_annotation[tmp_deleterious]:
                    tmp_deleterious=element_annotation_category
                else:
                    continue
        dicts[key][val_key]=tmp_deleterious


#returns variant for the gene and annoation(LOF1,LOF2...) using the dictionary
def get_key(gene, annot):
    list_key=[]
    for key, value in dicts.items():
        for value_key, value_value in value.items():
            if value_key == gene and value_value in annot:
                list_key.append(key)
    return list_key

#list_annot_run=list_annot_run2
#list_annot_run=list_annot_run2
dicts_skato={}
list_gene=[]
df_regenie=pd.read_csv(args.regenie_output, sep='\s+')

print("regenie summary statistics .tsv file from step9: \n")
print(df_regenie)

for row in df_regenie.itertuples():
    tmp_row=getattr(row, 'ID')
    tmp_gene=tmp_row.split('.')[0]
    list_gene.append(tmp_gene)

df_regenie['GENE']=list_gene
condition_mask1 = df_regenie['ID'].str.contains('Mask1')
condition_mask2 = df_regenie['ID'].str.contains('Mask2')
condition_mask3 = df_regenie['ID'].str.contains('Mask3')
condition_mask4 = df_regenie['ID'].str.contains('Mask4')
condition_af01 = df_regenie['ID'].str.contains('0.01')
print(list_annot_run == "list_annot_run1") 

if list_annot_run == "list_annot_run1":
    print('RUNNING run1')
    combined_condition_af01 = condition_mask1 & condition_af01
    df_regenie_Mask_af01 = df_regenie[combined_condition_af01]
    list_annot_run=list_annot_run1

elif list_annot_run == "list_annot_run2":
    print('RUNNING run2')
    combined_condition_af01 = condition_mask2 & condition_af01
    df_regenie_Mask_af01 = df_regenie[combined_condition_af01]
    list_annot_run=list_annot_run2

elif list_annot_run == "list_annot_run3":
    print('RUNNING run3')
    combined_condition_af01 = condition_mask3 & condition_af01
    df_regenie_Mask_af01 = df_regenie[combined_condition_af01]
    list_annot_run=list_annot_run3

elif list_annot_run == "list_annot_run4":
    print('RUNNING run4')
    combined_condition_af01 = condition_mask4 & condition_af01
    df_regenie_Mask_af01 = df_regenie[combined_condition_af01]
    list_annot_run=list_annot_run4

print("regenie summary statistics for selected Mask: \n")

print(df_regenie_Mask_af01)


#makes a new dictionary dicts_skato (gene:variant) by taking in the GENE and list_annot_run(LOF1..)
for row in df_regenie_Mask_af01.itertuples():
    tmp_id=getattr(row, 'GENE')
    #tmp_setID, tmp_Mask, tmp_afrq=tmp_id.rstrip().split('.',2)
    #print(tmp_setID) 
    #print(tmp_Mask)
    #print(tmp_afrq)
    list_gene=[]
    list_gene=get_key(tmp_id, list_annot_run)
    if tmp_id not in dicts_skato:
        dicts_skato[tmp_id]=list_gene
    else:
        print('extists! error')
print("done!")
print("dictionary for regenie with selected Mask: \n")
print(dicts_skato)


output_step2=pd.DataFrame()
for annot in list_annot_run:
    print("annot: ",annot)
    out_step1='list_'+annot+suffix
    #list_annot_directory='/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME/step2/'
    out1_step1=os.path.join(list_annot_directory, out_step1)
    print(out1_step1)
    df_run=pd.read_csv(out1_step1,sep='\t')
    df_run.columns=['SNP']
    df_check_run=pd.DataFrame()
    list_check_run=[]
    for row in df_run.itertuples():
        tmp_snp=getattr(row, 'SNP')
        tmp_chr, tmp_pos, tmp_ref, tmp_alt=tmp_snp.split(':')
        check_snp=tmp_snp
        if check_snp in dicts_GT_count_geneset[annot]:
            list_check_run.append(tmp_snp)
    df_check_run['SNP']=list_check_run
    summary_GT_case_control_run=pd.merge(summary_GT_case_control,df_check_run,on=['SNP'],how='inner')
    
    new_df=pd.DataFrame()
    list_gene=[]
    list_case_GT=[]
    list_control_GT=[]
    list_case_nonMiss=[]
    list_control_nonMiss=[]
    list_case_individual=[]
    list_control_individual=[]
    for key, value in dicts_skato.items():
        subset=summary_GT_case_control_run.loc[summary_GT_case_control_run["SNP"].isin(value)]
        if subset.empty == False:
            total=subset.sum()
            list_gene.append(key)
            list_control_GT.append(total[1])
            list_control_nonMiss.append(total[2])
            list_control_individual.append(len(total[3]))
            list_case_individual.append(len(total[4]))
            list_case_GT.append(total[5])
            list_case_nonMiss.append(total[6])
        elif subset.empty == True:
            list_gene.append(key)
            list_control_GT.append(pd.NA)
            list_control_nonMiss.append(pd.NA)
            list_control_individual.append(pd.NA)
            list_case_individual.append(pd.NA)
            list_case_GT.append(pd.NA)
            list_case_nonMiss.append(pd.NA)
    new_df['SNP']=list_gene
    new_df['CountGT_control_'+annot]=list_control_GT
    new_df['CountGT_case_'+annot]=list_case_GT
    new_df['nonMiss_control_'+annot]=list_control_nonMiss
    new_df['nonMiss_case_'+annot]=list_case_nonMiss
    new_df['#individual_allele_case_'+annot]=list_case_individual
    new_df['#individual_allele_control_'+annot]=list_control_individual
    if len(output_step2)==0:
        output_step2=new_df
    elif len(output_step2)>0:
        output_step2=pd.merge(output_step2, new_df, on=['SNP'], how='inner')
print("output_step2: \n")
print(output_step2)


prefix_control_GT='CountGT_control_'
prefix_case_GT='CountGT_case_'
prefix_control_nonMiss='nonMiss_control_'
prefix_case_nonMiss='nonMiss_case_'
output_step2_control_GT=output_step2.filter(like=prefix_control_GT, axis=1)
print(output_step2_control_GT)

output_step2_case_GT=output_step2.filter(like=prefix_case_GT, axis=1)
output_step2_control_nonMiss=output_step2.filter(like=prefix_control_nonMiss, axis=1)
output_step2_case_nonMiss=output_step2.filter(like=prefix_case_nonMiss, axis=1)
output_step2['CountGT_control_All']=output_step2_control_GT.sum(axis=1, skipna=True)
print(output_step2)

output_step2['CountGT_case_All']=output_step2_case_GT.sum(axis=1, skipna=True)
output_step2['CountGT_control_nonMiss']=output_step2_control_nonMiss.sum(axis=1, skipna=True)
output_step2['CountGT_case_nonMiss']=output_step2_case_nonMiss.sum(axis=1, skipna=True)
output_step2['Odds_Ratio']=((output_step2['CountGT_case_All']+1)*(output_step2['CountGT_control_nonMiss']-output_step2['CountGT_control_All']+1))/((output_step2['CountGT_control_All']+1)*(output_step2['CountGT_case_nonMiss']-output_step2['CountGT_case_All']+1))

array_chisquare=output_step2.loc[:,['CountGT_case_All','CountGT_case_nonMiss','CountGT_control_All','CountGT_control_nonMiss']]
array_chisquare['CountGT_case_All_other']=array_chisquare['CountGT_case_nonMiss']-array_chisquare['CountGT_case_All']
array_chisquare['CountGT_control_All_other']=array_chisquare['CountGT_control_nonMiss']-array_chisquare['CountGT_control_All']
contingency_table=array_chisquare.loc[:,['CountGT_case_All','CountGT_case_All_other','CountGT_control_All','CountGT_control_All_other']]

def chi_square_test(row):
    obs = [[row['CountGT_case_All'], row['CountGT_case_All_other']],
          [row['CountGT_control_All'], row['CountGT_control_All_other']]]
    chi2, p, dof, ex = chi2_contingency(obs)
    return p

p=contingency_table.apply(chi_square_test, axis=1)
output_step2['P_Chi-Square']=p
#Dir=/data/Segre_Lab/users/jlama/WES_new.ALL_050824/GeneBurden/FAME/step9_to_11/Run2_Count
Dir=os.path.join(args.output_path, args.output_name + ".tsv")
print(Dir)
output_step2.to_csv(Dir, sep='\t', index=None)
