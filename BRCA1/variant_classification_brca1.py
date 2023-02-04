import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from pytorch_tabnet.tab_model import TabNetClassifier
import torch
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, precision_recall_curve, roc_curve, roc_auc_score, matthews_corrcoef 
import seaborn as sns
from pandas import read_csv
from sklearn.impute import SimpleImputer
import numpy as np


#load exonic data
data = read_csv('pipeline.csv')
df = data.copy()

ID = df['ID']
df = df.drop('ID', 1)   

#split into input and output columns
X, y = df.iloc[:, :-1], df.iloc[:, -1]

#load the model
model2 = xgb.XGBClassifier()
model2.load_model("model_brca1.txt")
preds = model2.predict(X)

#save the predictions
df_preds = pd.DataFrame(preds)
df_preds.insert(0, 'ID', data['ID'])
header = ["ID", "Prediction"]
df_preds.columns = header
df_preds.to_csv("predictions_brca1.csv", index= False)
