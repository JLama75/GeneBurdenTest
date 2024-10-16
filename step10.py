import os
import sys
import numpy as np
import pandas as pd

df_run4_gb=pd.read_csv('./step8_to_10/Regenie_Run4.tsv', sep='\t')
df_run4_gb=df_run4_gb[['GENE','P']]
df_run4_gb.columns=['GENE','P-value.Run4']
#print(df_run4_gb)

df_run1_gb=pd.read_csv('./step8_to_10/Regenie_Run1.tsv', sep='\t')
df_run1_step2=pd.read_csv('./step8_to_10/Run1_Count.tsv', sep='\t')
df_run1=pd.merge(df_run1_gb, df_run1_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run1=df_run1.drop(['SNP'], axis=1)
df_run1=pd.merge(df_run1,df_run4_gb, on=['GENE'], how='inner')

print(df_run1)
df_run1.to_csv('./step8_to_10/Run1_merged_Count.tsv', index=None, sep='\t')

df_run2_gb=pd.read_csv('./step8_to_10/Regenie_Run2.tsv', sep='\t')
df_run2_step2=pd.read_csv('./step8_to_10/Run2_Count.tsv', sep='\t')
df_run2=pd.merge(df_run2_gb, df_run2_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run2=df_run2.drop(['SNP'], axis=1)
df_run2=pd.merge(df_run2,df_run4_gb, on=['GENE'], how='inner')
print(df_run2)
df_run2.to_csv('./step8_to_10/Run2_merged_Count.tsv', index=None, sep='\t')

df_run3_gb=pd.read_csv('./step8_to_10/Regenie_Run3.tsv', sep='\t')
df_run3_step2=pd.read_csv('./step8_to_10/Run3_Count.tsv', sep='\t')
df_run3=pd.merge(df_run3_gb, df_run3_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run3=df_run3.drop(['SNP'], axis=1)
df_run3=pd.merge(df_run3,df_run4_gb, on=['GENE'], how='inner')
print(df_run3)
df_run3.to_csv('./step8_to_10/Run3_merged_Count.tsv', index=None, sep='\t')

df_run4_gb=pd.read_csv('./step8_to_10/Regenie_Run4.tsv', sep='\t')
df_run4_step2=pd.read_csv('./step8_to_10/Run4_Count.tsv', sep='\t')
df_run4=pd.merge(df_run4_gb, df_run4_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run4=df_run4.drop(['SNP'], axis=1)
print(df_run4)
df_run4.to_csv('./step8_to_10/Run4_merged_Count.tsv', index=None, sep='\t')
