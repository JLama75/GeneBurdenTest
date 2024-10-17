import os
import sys
import numpy as np
import pandas as pd
import argparse
from statsmodels.stats.multitest import multipletests
import xlsxwriter

#Set up argument parser
parser = argparse.ArgumentParser(description="Final data preparation for each Runs")

# Add arguments
parser.add_argument('--Dir', required=True, help="Path to the input foders")

# Parse the arguments
args = parser.parse_args()

Dir=args.Dir #./step8_to_10/

Count1=os.path.join(Dir, "Run1_Count.tsv")
Count2=os.path.join(Dir, "Run2_Count.tsv")
Count3=os.path.join(Dir, "Run3_Count.tsv")
Count4=os.path.join(Dir, "Run4_Count.tsv")

Regenie1=os.path.join(Dir, "Regenie_Run1.tsv")
Regenie2=os.path.join(Dir, "Regenie_Run2.tsv")
Regenie3=os.path.join(Dir, "Regenie_Run3.tsv")
Regenie4=os.path.join(Dir, "Regenie_Run4.tsv")


Out1=os.path.join(Dir, "Run1_merged_Count.tsv")
Out2=os.path.join(Dir, "Run2_merged_Count.tsv")
Out3=os.path.join(Dir, "Run3_merged_Count.tsv")
Out4=os.path.join(Dir, "Run4_merged_Count.tsv")

excelOut=os.path.join(Dir, "SummaryTable_allRuns.xlsx")

def calculate_FDR(df, p_value_column, new_column='FDR_BH_p_value'):
    """
    Applies the Benjamini-Hochberg FDR correction to a DataFrame's P column.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the P-values.
    p_value_column (str): The name of the column containing the raw P-values.
    new_column (str): The name of the new column where FDR-corrected P-values will be stored.
                      Default is 'FDR_BH_p_value'.

    Returns:
    pd.DataFrame: The DataFrame with an additional column for FDR-corrected P-values.
    """
    # Perform Benjamini-Hochberg correction
    df[new_column] = multipletests(df[p_value_column], method='fdr_bh')[1]
    
    # Return the updated DataFrame
    return df


df_run4_gb=pd.read_csv(Regenie4, sep='\t')
df_run4_gb=df_run4_gb[['GENE','P']]
df_run4_gb.columns=['GENE','P-value.Run4']
#print(df_run4_gb)

df_run1_gb=pd.read_csv(Regenie1, sep='\t')
df_run1_step2=pd.read_csv(Count1, sep='\t')
df_run1=pd.merge(df_run1_gb, df_run1_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run1=df_run1.drop(['SNP'], axis=1)
df_run1=pd.merge(df_run1,df_run4_gb, on=['GENE'], how='left')
df_run1=calculate_FDR(df_run1, 'P')
df_run1.fillna('NA', inplace=True)

print(df_run1)
print(f"Number of rows for run1: {df_run1.shape[0]}")
print(f"Number of cols for run1: {df_run1.shape[1]}")

df_run1.to_csv(Out1, index=None, sep='\t')

df_run2_gb=pd.read_csv(Regenie2, sep='\t')
df_run2_step2=pd.read_csv(Count2, sep='\t')
df_run2=pd.merge(df_run2_gb, df_run2_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run2=df_run2.drop(['SNP'], axis=1)
df_run2=pd.merge(df_run2,df_run4_gb, on=['GENE'], how='left')
df_run2=calculate_FDR(df_run2, 'P')
df_run2.fillna('NA', inplace=True)

print(df_run2)
print(f"Number of rows for run2: {df_run2.shape[0]}")
print(f"Number of cols for run2: {df_run2.shape[1]}")
df_run2.to_csv(Out2, index=None, sep='\t')

df_run3_gb=pd.read_csv(Regenie3, sep='\t')
df_run3_step2=pd.read_csv(Count3, sep='\t')
df_run3=pd.merge(df_run3_gb, df_run3_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run3=df_run3.drop(['SNP'], axis=1)
df_run3=pd.merge(df_run3,df_run4_gb, on=['GENE'], how='left')
df_run3=calculate_FDR(df_run3, 'P')
df_run3.fillna('NA', inplace=True)

print(df_run3)
print(f"Number of rows for run3: {df_run3.shape[0]}")
print(f"Number of cols for run3: {df_run3.shape[1]}")
df_run3.to_csv(Out3, index=None, sep='\t')

df_run4_gb=pd.read_csv(Regenie4, sep='\t')
df_run4_step2=pd.read_csv(Count4, sep='\t')
df_run4=pd.merge(df_run4_gb, df_run4_step2, left_on=['GENE'], right_on=['SNP'], how='inner')
df_run4=df_run4.drop(['SNP'], axis=1)
df_run4=calculate_FDR(df_run4, 'P')
df_run4.fillna('NA', inplace=True)

print(df_run4)
print(f"Number of rows for run4: {df_run4.shape[0]}")
print(f"Number of cols for run4: {df_run4.shape[1]}")
df_run4.to_csv(Out4, index=None, sep='\t')

with pd.ExcelWriter(excelOut, engine='xlsxwriter') as writer:
    df_run1_gb.to_excel(writer, sheet_name='Run1_merged_Count', index=False)
    df_run2_gb.to_excel(writer, sheet_name='Run2_merged_Count', index=False)
    df_run3_gb.to_excel(writer, sheet_name='Run3_merged_Count', index=False)
    df_run4_gb.to_excel(writer, sheet_name='Run4_merged_Count', index=False)
    
    
