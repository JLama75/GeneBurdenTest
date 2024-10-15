import os
import sys
import pandas as pd
import numpy as np
import json
import argparse

f = open('./step2_to_4/Database_Gene_Variant_Annot.json')
data = json.load(f)

#define four list
df_annot=pd.DataFrame()

list_symbol=[]
list_snp=[]
list_annot=[]
#variant: i
#gene: gene
#annotation: data[i][gene]

# Iterating through the json
for i in data:
    for gene in data[i]:
        list_annot.append(data[i][gene])
        list_symbol.append(gene)
        list_snp.append(i)

df_annot['Symbol']=list_symbol
df_annot['SNP']=list_snp
df_annot['Annotation']=list_annot

# Closing file
f.close()

print(df_annot)
df_annot=df_annot.loc[:,['SNP','Symbol','Annotation']]
df_annot.to_csv('./step2_to_4/study_annot.tsv', index=None, sep='\t', header=None)
