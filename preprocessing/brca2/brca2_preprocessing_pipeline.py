from doctest import FAIL_FAST
import os
import subprocess

os.system("echo Hello from the other side!")


#create input file for pfam:
print("Annotating pfam")

#command will create bed file with ID as 4th column to enable bedtools intersect 
COMMAND = '''awk -F"\t" 'OFS="\t" {print "chr"$158, $159, $4, $1}' brca2_multianno_processed | grep -v "ID" | grep -v "START" > mydata.bed'''
subprocess.call(COMMAND, shell=True)
print("Bedfile created")

#Next, perform bedtools intersect to create input file for pfam to merge with brca2_multianno_processed 
COMMAND = '''bedtools intersect -a mydata.bed -b ../../pfam_all.bed -wa | awk -F"\t" 'OFS="\t" {print $4, "1"}' - | sort -u - | awk 'BEGIN{print "ID\tPfam_imp_domain"}1' - > op_pfam'''
subprocess.call(COMMAND, shell=True)
print("Intersect file created\n")


###READ ANNOVAR O/P AND PREPROCESS THE FILE###
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from pandas import read_csv
import numpy as np
import re

##preprocessing
#read the annovar processed o/p tsv file and convert it into a csv
data1 = pd.read_table('brca2_multianno_processed', sep='\t')
df = data1.copy()
print("The input file shape is: ")
print(df.shape)

#rename Func.refGene to Function and add Outcome column
df.rename(columns = {'ExonicFunc.refGene':'Function'}, inplace = True)
print("Renamed ExonicFunc to Function")

#keep only exonic variants:
df = df[df['Func.refGene'] == "exonic"]
print("All non-exonic variants removed ..")
print(df.shape)
print("Removed non exonic variants")

#label encoding the ExonicFunc.refGene column:
cleanup_nums = {"Function" : {"downstream" : 1, "nonsynonymous SNV" : 2, "frameshift insertion" : 3, "frameshift deletion" : 4, "stopgain" : 5, "nonframeshift deletion" : 6, "frameshift substitution" : 7, "startloss" : 8, "synonymous SNV" : 9, "nonframeshift substitution" : 10, "stoploss" : 11, "nonframeshift insertion" : 12, "intergenic" : 13, "intronic" : 14, "splicing" : 15, "UTR3" : 16, "UTR5" : 17, "upstream" : 18, "ncRNA_intronic" : 19}}
df = df.replace(cleanup_nums)
print('Labels encoded ..')
print(df.shape)


#merge with external result files:
outcome = read_csv("op_outcome", sep="\t")
df = pd.merge(df, outcome, on='ID', how='left')
print('Outcome successfully added ..')
print(df.shape)
loftee = read_csv("op_lof", sep = "\t")
df = df.merge(loftee, on='ID', how = 'left')
print('LoF successfully added ..')
print(df.shape)
pfam = read_csv("op_pfam", sep ="\t")
df = df.merge(pfam, on='ID', how = 'left')
print('Pfam domains added ..')
print(df.shape)


#extract Start_pro
df['Start_pro'] = df['AAChange.refGene'].str.extract('p\.[A-Z](\d+)[\w\.*|\W\.*,]')
print("Protein start positions extracted and added ..")
print(df.shape)
   
#rearrange column names according to final order:
df2 = read_csv('../../final_colnames.csv')
data2 = df2.columns.tolist()
column_list = data2

#shuffle columns based on file that trained our model
shuffled_new = df[column_list]
shuffled_new = df.reindex(columns=column_list)
print('Columns reindexed ..')
print(shuffled_new.shape)

# #replace dots but not decimals with nan
shuffled_new = shuffled_new.replace('(?<!.)\.(?!.)', np.NaN, regex=True)
print("Replaced dots with NaN ..")
#print(shuffled_new[:1])
print(shuffled_new.shape)

#replace NaN with 0 in AF columns
list_of_AFs = list(shuffled_new.loc[0:0, 'LoF_HC_Canonical':'SAS.sites.2015_08'])
shuffled_new[list_of_AFs] = shuffled_new[list_of_AFs].replace({np.nan:0, np.nan:0})
print("Replaced NaN with 0 for AFs ..")
print(shuffled_new.shape)

#save file
pd.DataFrame(shuffled_new).to_csv('../../BRCA2/pipeline.csv', index=False)


