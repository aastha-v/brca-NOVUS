import os
import subprocess
import xgboost as xgb
import pandas as pd
from pandas import read_csv
import numpy as np
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, precision_recall_curve, roc_curve, roc_auc_score, matthews_corrcoef 
import seaborn as sns
from numpy import nan
import sys

os.system("echo Hello from the other side!")

#create input file for pfam:
print("Annotating pfam")

#command will create bed file with ID as 4th column to enable bedtools intersect 
COMMAND = '''awk -F"\t" 'OFS="\t" {print "chr"$158, $159, $4, $1}' multianno_processed | grep -v "ID" | grep -v "START" > mydata.bed'''
subprocess.call(COMMAND, shell=True)
print("Bedfile created")

#Next, perform bedtools intersect to create input file for pfam to merge with multianno_all_processed in python
COMMAND = '''bedtools intersect -a mydata.bed -b ../scripts/pfam_all.bed -wa | awk -F"\t" 'OFS="\t" {print $4, "1"}' - | sort -u - | awk 'BEGIN{print "ID\tPfam_imp_domain"}1' - > op_pfam'''
subprocess.call(COMMAND, shell=True)
print("Intersect file created\n")

##preprocessing
#read the annovar processed o/p tsv file and convert it into a csv
data1 = pd.read_table('multianno_processed', sep='\t')
df = data1.copy()
print("The input file shape is: ")
print(df.shape)

#rename Func.refGene to Function and add Outcome column
df.rename(columns = {'ExonicFunc.refGene':'Function'}, inplace = True)
print("Renamed ExonicFunc to Function")

#keep only exonic variants:
df = df[df['Func.refGene'] == "exonic"]
print("All non-exonic variants removed ..")
print("Removed non exonic variants")

#label encoding the ExonicFunc.refGene column:
cleanup_nums = {"Function" : {"downstream" : 1, "nonsynonymous SNV" : 2, "frameshift insertion" : 3, "frameshift deletion" : 4, "stopgain" : 5, "nonframeshift deletion" : 6, "frameshift substitution" : 7, "startloss" : 8, "synonymous SNV" : 9, "nonframeshift substitution" : 10, "stoploss" : 11, "nonframeshift insertion" : 12, "intergenic" : 13, "intronic" : 14, "splicing" : 15, "UTR3" : 16, "UTR5" : 17, "upstream" : 18, "ncRNA_intronic" : 19}}
df = df.replace(cleanup_nums)
print('Labels encoded ..')

#merge with external result files:
outcome = read_csv("../input_folder/op_outcome_train_brca1", sep="\t")
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
df['Start_pro'] = df['AAChange.refGene'].str.extract('NM_007294:(.*),')
df['Start_pro'] = df['Start_pro'].str.extract('p\.[A-Z](\d+)[\w\.*|\W\.*,]')
print("Protein start positions extracted and added ..")
print(df.shape)
   
#rearrange column names according to final order:
df2 = read_csv('../scripts/final_colnames.csv')
data2 = df2.columns.tolist()
column_list = data2

#shuffle columns based on file that trained our model
shuffled_new = df[column_list]
shuffled_new = df.reindex(columns=column_list)
print('Columns reindexed ..')
print(shuffled_new.shape)

#replace dots but not decimals with nan
shuffled_new = shuffled_new.replace('(?<!.)\.(?!.)', np.NaN, regex=True)
print("Replaced dots with NaN ..")

#replace NaN with 0 in AF columns
list_of_AFs = list(shuffled_new.loc[0:0, 'LoF_HC_Canonical':'SAS.sites.2015_08'])
shuffled_new[list_of_AFs] = shuffled_new[list_of_AFs].replace({np.nan:0, np.nan:0})
print("Replaced NaN with 0 for AFs ..")

#save file
pd.DataFrame(shuffled_new).to_csv('pipeline.csv', index=False)

#model creation and use
var = sys.argv[1]

if var == "create":
    #model creation
    #load all data
    data = read_csv('pipeline.csv')
    df3 = data.copy()
    df3 = df.drop('ID', 1)
    
    #subset the data into x and y
    X, y = df.iloc[:,:-1],df.iloc[:,-1]
    
    
    #split data into test and train
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    
    model = xgb.XGBClassifier(scale_pos_weight=3.888099467, base_score=0.5, booster='gbtree', callbacks=None,
                  colsample_bylevel=1, colsample_bynode=1, colsample_bytree=0.3,
                  early_stopping_rounds=None, enable_categorical=False,
                  eval_metric=None, gamma=0.4, gpu_id=-1, grow_policy='depthwise',
                  importance_type=None, interaction_constraints='',
                  learning_rate=0.1, max_bin=256, max_cat_to_onehot=4,
                  max_delta_step=0, max_depth=4, max_leaves=0, min_child_weight=1,
                  monotone_constraints='()', n_estimators=150,
                  n_jobs=0, num_parallel_tree=1, predictor='auto', random_state=0,
                  reg_alpha=0, reg_lambda=1, missing=nan, seed=123)
    
    #fit the model and make predictions
    model.fit(X_train,y_train)
    preds = model.predict(X_test)
    
    #save the model
    model.save_model("../model_creation/selfmodel_brca1.txt")
    
    
    ###CHECK MODEL PERFORMANCE
    # plot auc
    fig3 = plt.figure()
    y_pred_proba = model.predict_proba(X_test)[:,1]
    fpr, tpr, threshold = roc_curve(y_test,  y_pred_proba)
    
    #Compute Area Under the Receiver Operating Characteristic Curve (ROC AUC) from prediction scores
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    plt.plot(fpr,tpr)
    plt.title('ROC Curve')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig("../model_creation/brca1/auc.png")
    
    
    # plot precision recall curve
    fig3_2 = plt.figure()
    precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
    plt.plot(recall, precision)
    plt.title('Precision-Recall Curve')
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    plt.savefig("../model_creation/brca1/prcurve.png")
    
     
    # determine best accuracy for test set
    preds = model.predict(X_test)
    test_acc = accuracy_score(preds, y_test)
    
    #determine matthews correlation coefficient
    MCC = matthews_corrcoef(y_test, preds)
    
    print(f"BEST ACCURACY SCORE ON TEST SET : {test_acc}")
    print("AUC ROC : ", roc_auc)
    print("MCC : ", MCC)
    
    print(classification_report(y_test, preds))
    report = classification_report(y_test, preds, output_dict=True)
    report_df = pd.DataFrame(report).transpose()
    report_df.to_csv("../model_creation/brca1/classification_report.csv")
    
    #plot confusion matrix
    fig5 = plt.figure()
    sns.heatmap(confusion_matrix(y_test, preds), cmap = 'Blues', annot = True, fmt = 'd', linewidths = 5, cbar = False, annot_kws = {'fontsize': 15},
                yticklabels = ['Benign', 'Pathogenic'], xticklabels = ['Predicted Benign', 'Predicted Pathogenic'])
    #plt.yticks(rotation = 0)
    plt.title('Confusion Matrix')
    plt.savefig("../model_creation/brca1/confusion_matrix.png")
    
    
    #Visualize Boosting Trees and Feature Importance using plot_tree() function
    #boosting trees
    figure6 = plt.figure()
    xgb.plot_tree(model,num_trees=0)
    plt.rcParams['figure.figsize'] = [100, 200]
    plt.savefig("../model_creation/brca1/boosted_trees.png")
    
    #feature importance
    figure7 = plt.figure()
    xgb.plot_importance(model)
    plt.rcParams['figure.figsize'] = [7, 7]
    plt.margins()
    plt.savefig("../model_creation/brca1/feature_importance.png", bbox_inches='tight')
    
else:
    #load exonic data
    data = read_csv('pipeline2.csv')
    df4 = data.copy()
    
    ID = df4['ID']
    df4 = df4.drop('ID', 1)   
    
    #split into input and output columns
    X, y = df.iloc[:, :-1], df.iloc[:, -1]
    
    
    #load the model
    model2 = xgb.XGBClassifier()
    model2.load_model("../model_creation/selfmodel_brca1.txt")
    #preds = model2.predict(X)
    
    #save the predictions
    #df_preds = pd.DataFrame(preds)
    #df_preds.insert(0, 'ID', data['ID'])
    #header = ["ID", "Prediction"]
    #df_preds.columns = header
    #df_preds.to_csv("../model_creation/brca1/selfpredictions_brca1.csv", index= False)
    print("kthnbye")
