    
import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from pandas import read_csv
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, precision_recall_curve, roc_curve, roc_auc_score, matthews_corrcoef 
import seaborn as sns
from numpy import nan


#load all data
data = read_csv('brca2_pipeline.csv')
df = data.copy()
df = df.drop('ID', 1)

#subset the daya into x and y
X, y = df.iloc[:,:-1],df.iloc[:,-1]

#split data into test and train
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=123)

model = xgb.XGBClassifier(scale_pos_weight=3.046349942, base_score=0.5, booster='gbtree', callbacks=None,
              colsample_bylevel=1, colsample_bynode=1, colsample_bytree=0.3,
              early_stopping_rounds=None, enable_categorical=False,
              eval_metric=None, gamma=0.3, gpu_id=-1, grow_policy='depthwise',
              importance_type=None, interaction_constraints='',
              learning_rate=0.2, max_bin=256, max_cat_to_onehot=4,
              max_delta_step=0, max_depth=15, max_leaves=0, min_child_weight=5,
              monotone_constraints='()', n_estimators=50, n_jobs=0,
              num_parallel_tree=1, predictor='auto', random_state=0,
              reg_alpha=0, reg_lambda=1, missing=nan, seed=124)

#fit the model and make predictions
model.fit(X_train,y_train)
preds = model.predict(X_test)

#save the model
model.save_model("model_brca2.txt")


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
plt.savefig("brca2/auc.png", dpi=600)


# plot precision recall curve
fig3_2 = plt.figure()
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
plt.plot(recall, precision)
plt.title('Precision-Recall Curve')
plt.ylabel('Precision')
plt.xlabel('Recall')
plt.savefig("brca2/prcurve.png", dpi=600)


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
report_df.to_csv("brca2/classification_report.csv")

#plot confusion matrix
fig5 = plt.figure()
sns.heatmap(confusion_matrix(y_test, preds), cmap = 'Blues', annot = True, fmt = 'd', linewidths = 5, cbar = False, annot_kws = {'fontsize': 15},
            yticklabels = ['Benign', 'Pathogenic'], xticklabels = ['Predicted Benign', 'Predicted Pathogenic'])
plt.title('Confusion Matrix')
plt.savefig("brca2/confusion_matrix.png", dpi=600)


#Visualize Boosting Trees and Feature Importance using plot_tree() function
#boosting trees
figure6 = plt.figure()
xgb.plot_tree(model,num_trees=0)
plt.rcParams['figure.figsize'] = [100, 200]
plt.savefig("brca2/boosted_trees.png", dpi=600)

#feature importance
figure7 = plt.figure()
xgb.plot_importance(model)
plt.rcParams['figure.figsize'] = [7, 7]
plt.margins()
plt.savefig("brca2/feature_importance.png", bbox_inches='tight', dpi=600)

