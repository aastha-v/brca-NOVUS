import xgboost as xgb
import pandas as pd
from pandas import read_csv


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
