import os
import sys
import pandas as pd
import numpy as np
import json
import argparse

df=pd.read_csv('./step2_to_4/study_annot.tsv', sep='\t', header=None)
df.columns=['SNP','Symbol','Annotation']
print(df)
df_combined = df.groupby('Symbol')['SNP'].agg(lambda x: ','.join(x)).reset_index()

list_chr=[]
list_pos=[]
for row in df_combined.itertuples():
    tmp_SNP=getattr(row, 'SNP')
    tmp_SNP_1st=tmp_SNP.split(',',1)[0]
    tmp_chr, tmp_pos, tmp_ref, tmp_alt=tmp_SNP_1st.split(':')
    list_chr.append(tmp_chr)
    list_pos.append(tmp_pos)
df_combined['CHR']=list_chr
df_combined['POS']=list_pos
print(df_combined)

df_combined=df_combined.loc[:,['Symbol','CHR','POS','SNP']]
#df_combined.to_csv('./step2_to_4/study_set.tsv', index=None, sep='\t', header=False)

# Function to check SNP matches CHR
def check_chr_match(row):
    chr_val = str(row['CHR'])  # Convert CHR to string for comparison
    snps = row['SNP'].split(',')  # Split the SNP string by commas

    for snp in snps:
        snp_chr = snp.split(':')[0]  # Extract the chromosome part from SNP
        if snp_chr != chr_val:  # Check if SNP chromosome matches CHR
            return False  # If any SNP doesn't match, return False

    return True  # Return True if all SNPs match

# Apply the function and split DataFrame into matched and non-matched
df_matched = df_combined[df_combined.apply(check_chr_match, axis=1)]
df_not_matched = df_combined[~df_combined.apply(check_chr_match, axis=1)]
print(df_not_matched)

df_matched.to_csv('./step2_to_4/study_set.tsv', index=None, sep='\t', header=False)
df_not_matched.to_csv('./step2_to_4/removedIDs.studySet.tsv', index=None, sep='\t', header=False)

# Select the "Symbol" column
#df_not_matched = df_not_matched['Symbol']
# Save the column to a .tsv file without the column header
#df_not_matched.to_csv('./step2_to_4/exclude.study_set.tsv', sep='\t', header=False, index=False)

f = open('./step2_to_4/Database_Gene_Variant_Annot.json')
data = json.load(f)

#define four list
df_LoF1=pd.DataFrame()
df_LoF2=pd.DataFrame()
df_missense=pd.DataFrame()
df_moderate=pd.DataFrame()
df_modifier=pd.DataFrame()
df_low=pd.DataFrame()
df_synonymous=pd.DataFrame()

list_symbol_LoF1=[]
list_snp_LoF1=[]
list_symbol_LoF2=[]
list_snp_LoF2=[]
list_symbol_missense=[]
list_snp_missense=[]
list_symbol_moderate=[]
list_snp_moderate=[]
list_symbol_modifier=[]
list_snp_modifier=[]
list_symbol_low=[]
list_snp_low=[]
list_symbol_synonymous=[]
list_snp_synonymous=[]

#variant: i
#gene: gene
#annotation: data[i][gene]

# Iterating through the json
for i in data:
    for gene in data[i]:
        if data[i][gene] == 'LoF1':
            list_symbol_LoF1.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_LoF1.append(i)
        if data[i][gene] == 'LoF2':
            list_symbol_LoF2.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_LoF2.append(i)
        if data[i][gene] == 'missense':
            list_symbol_missense.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_missense.append(i)
        if data[i][gene] == 'moderate':
            list_symbol_moderate.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_moderate.append(i)
        if data[i][gene] == 'modifier':
            list_symbol_modifier.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_modifier.append(i)
        if data[i][gene] == 'low':
            list_symbol_low.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_low.append(i)
        if data[i][gene] == 'synonymous':
            list_symbol_synonymous.append(gene)
#             tmp_chr, tmp_pos, tmp_ref, tmp_alt=i.rstrip().split(':')
#             new_snp=':'.join([tmp_chr[3:], tmp_pos])
            list_snp_synonymous.append(i)
df_LoF1['Symbol']=list_symbol_LoF1
df_LoF1['SNP']=list_snp_LoF1

df_LoF2['Symbol']=list_symbol_LoF2
df_LoF2['SNP']=list_snp_LoF2

df_missense['Symbol']=list_symbol_missense
df_missense['SNP']=list_snp_missense

df_moderate['Symbol']=list_symbol_moderate
df_moderate['SNP']=list_snp_moderate

df_modifier['Symbol']=list_symbol_modifier
df_modifier['SNP']=list_snp_modifier

df_low['Symbol']=list_symbol_low
df_low['SNP']=list_snp_low

df_synonymous['Symbol']=list_symbol_synonymous
df_synonymous['SNP']=list_snp_synonymous

# Closing file
f.close()
print("print LoF1 dataframe: \n")
print(df_LoF1)
df_LoF1=df_LoF1.drop_duplicates()
df_LoF2=df_LoF2.drop_duplicates()
df_missense=df_missense.drop_duplicates()
df_moderate=df_moderate.drop_duplicates()
df_modifier=df_modifier.drop_duplicates()
df_low=df_low.drop_duplicates()
df_synonymous=df_synonymous.drop_duplicates()

df_LoF1.to_csv('./step2_to_4/study_LoF1.SetID', sep='\t', header=False, index=None)
df_LoF2.to_csv('./step2_to_4/study_LoF2.SetID', sep='\t', header=False, index=None)
df_missense.to_csv('./step2_to_4/study_missense.SetID', sep='\t', header=False, index=None)
df_moderate.to_csv('./step2_to_4/study_moderate.SetID', sep='\t', header=False, index=None)
df_modifier.to_csv('./step2_to_4/study_modifier.SetID', sep='\t', header=False, index=None)
df_low.to_csv('./step2_to_4/study_low.SetID', sep='\t', header=False, index=None)
df_synonymous.to_csv('./step2_to_4/study_synonymous.SetID', sep='\t', header=False, index=None)
