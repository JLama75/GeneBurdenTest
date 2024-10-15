import os
import sys
import json
import numpy as np
import argparse
import pandas as pd
import re
from collections import defaultdict
import json
import subprocess
from collections import Counter

file_vcf = sys.argv[1]
DIR = './step2_to_4/'
print(f"Directory path: {DIR}")
print(f"File name: {file_vcf}")

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

def Merge(dict1, dict2):
    for i in dict2.keys():
        if i in dict1:
            if dicts_annotation[dict1[i]] <= dicts_annotation[dict2[i]]:
                dict1[i] = dict1[i]
            elif dicts_annotation[dict1[i]] > dicts_annotation[dict2[i]]:
                dict1[i] = dict2[i]
        elif i not in dict1:
            dict1[i] = dict2[i]
    return dict1

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

#create a dictionary with value as a list
dicts=nested_dict(2,list)
#with open(file_vcf_out) as vcf_in, open(file_vcf_out1, 'w') as vcf_out:
count=0
with open(file_vcf) as vcf_in:
    for line in vcf_in:
        if line[0]=='#':
           continue
        tmp_Variation,tmp_Location,tmp_Allele,tmp_Gene,tmp_Feature,tmp_Feature_type,tmp_Consequence,tmp_cDNA_position,tmp_CDS_position,tmp_Protein_position,tmp_Amino_acids,tmp_Codons,tmp_Existing_variation,tmp_Extra = line.rstrip().split('\t')
        tmp_annot=tmp_Consequence
        count=count+1
        if 'SYMBOL=' in tmp_Extra:
           result=re.search('SYMBOL=(.*);SYMBOL_SOURCE', tmp_Extra)
           tmp_gene=result.group(1)
        elif 'SYMBOL=' not in tmp_Extra:
           continue
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
print("count the lines in vcf: ",count)
variant_count = len(dicts)
print("count the no. of variants in dictionary dicts: ", variant_count)

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
variant_count = len(dicts)
print("count the no. of variants in dictionary dicts 2: ", variant_count)

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
        
variant_count = len(dicts)
print("count the no. of variants in dictionary dicts: ", variant_count)

#save dictionary to JSON file
#DIR = os.path.join(directory_path, "step2_to_4")
JSON = os.path.join(DIR, "Database_Gene_Variant_Annot.json")
print(JSON)
with open(JSON,'w') as json_file:
   json.dump(dicts, json_file)

list_LoF1=[]
list_LoF2=[]
list_missense=[]
list_moderate=[]
list_modifier=[]
list_low=[]
list_synonymous=[]
for key, val in dicts.items():
    for val_key, val_val in val.items():
        if val_val=='LoF1':
            list_LoF1.append(key)
        elif val_val=='LoF2':
            list_LoF2.append(key)
        elif val_val=='missense':
            list_missense.append(key)
        elif val_val=='moderate':
            list_moderate.append(key)
        elif val_val=='modifier':
            list_modifier.append(key)
        elif val_val=='low':
            list_low.append(key)
        elif val_val=='synonymous':
            list_synonymous.append(key)

df_LoF1=pd.DataFrame()
df_LoF2=pd.DataFrame()
df_missense=pd.DataFrame()
df_moderate=pd.DataFrame()
df_modifier=pd.DataFrame()
df_low=pd.DataFrame()
df_synonymous=pd.DataFrame()

df_LoF1['variant']=list_LoF1
df_LoF2['variant']=list_LoF2
df_missense['variant']=list_missense
df_moderate['variant']=list_moderate
df_modifier['variant']=list_modifier
df_low['variant']=list_low
df_synonymous['variant']=list_synonymous

LoF1 = os.path.join(DIR, "list_LoF1.tsv")
LoF2 = os.path.join(DIR, "list_LoF2.tsv")
missense = os.path.join(DIR, "list_missense.tsv")
moderate = os.path.join(DIR, "list_moderate.tsv")
modifier = os.path.join(DIR, "list_modifier.tsv")
low = os.path.join(DIR, "list_low.tsv")
synonymous = os.path.join(DIR, "list_synonymous.tsv")

df_LoF1.to_csv(LoF1, index=None, sep='\t')
df_LoF2.to_csv(LoF2, index=None, sep='\t')
df_missense.to_csv(missense, index=None, sep='\t')
df_moderate.to_csv(moderate, index=None, sep='\t')
df_modifier.to_csv(modifier, index=None, sep='\t')
df_low.to_csv(low, index=None, sep='\t')
df_synonymous.to_csv(synonymous, index=None, sep='\t')
