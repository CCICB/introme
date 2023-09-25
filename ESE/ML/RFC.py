# A Random Forest Classifier determining Splice Altering or Normal
# Import modules
import pandas as pd
CLASSIFICATION = True
if CLASSIFICATION:
    from sklearn.ensemble import RandomForestClassifier
else:
    from sklearn.ensemble import RandomForestRegressor

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import numpy as np
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import average_precision_score, precision_recall_curve, PrecisionRecallDisplay
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error, f1_score
import matplotlib.pyplot as plt
from sklearn.feature_selection import RFECV
import sys
import pickle

# #### LOAD DATA ####
data = sys.argv[1] # vcfAnno/data/output/neg.tsv
output = open(sys.argv[2], 'w+')
# pos_file = sys.argv[2] # vcfAnno/data/output/pos.tsv

# # Train data 
df = pd.read_csv(data, sep='\t', na_values=["."])
if CLASSIFICATION:
    # Drop "Low-freq" rows
    df = df[df['classification'] != 'Low-freq']
    df = pd.get_dummies(df, columns=['classification', 'location'], dtype=int)

    # Select features and target
    features = df.drop(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'gene_id', 'strand', 'MFASS_delta_index',
                                'classification_Normal', 'location_Intronic', 'classification_Splice-altering'])
    target = df['classification_Splice-altering']

else:
    df = pd.get_dummies(df, columns=['location'], dtype=int)
    features = df.drop(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'gene_id', 'strand', 'MFASS_delta_index',
                                'classification', 'location_Intronic'])
    target = df['MFASS_delta_index']

print(features)
print(target)
features.to_csv(output, encoding='utf-8', index=False, sep='\t')

X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.15, random_state=42)

# Initializing and training the RandomForestRegressor
if CLASSIFICATION:
    regressor = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=None, min_samples_split=3, min_samples_leaf=2, criterion='log_loss', class_weight='balanced_subsample', n_jobs=-1)
    regressor.fit(X_train, y_train)
else:
    regressor = RandomForestRegressor(n_estimators=400, random_state=42, max_depth=171, max_features='sqrt', criterion='squared_error', n_jobs=-1)
    regressor.fit(X_train, y_train)

    y_pred2 = regressor.predict(X_train)
    mae = mean_absolute_error(y_train, y_pred2)
    mse = mean_squared_error(y_train, y_pred2)
    rmse = mean_squared_error(y_train, y_pred2, squared=False)
    r2 = r2_score(y_train, y_pred2)


    print(f"{mae=}, {mse=}, {r2=}, {rmse=}")

    # Predicting and evaluating the model
    y_pred = regressor.predict(X_test)
    mae = mean_absolute_error(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)


    print(f"{mae=}, {mse=}, {r2=}, {rmse=}")
    exit(0)

# Predict probabilities
y_scores = regressor.predict_proba(X_test)[:, 1]

# Compute Precision, Recall
precision, recall, thresholds = precision_recall_curve(y_test, y_scores)

# Compute Average Precision (AP)
ap = average_precision_score(y_test, y_scores)

# Plot PR curve
plt.figure(figsize=(10, 7))
plt.plot(recall, precision, label=f'AP (Average Precision) = {ap:.2f}')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()
plt.grid(True)
# plt.show()
plt.savefig("RCRUNCH_log_e200_d170.png")

# Compute F1 score (at a specific threshold, e.g., 0.5 for binary classification)
y_pred = regressor.predict(X_test)
f1 = f1_score(y_test, y_pred)
print(f'F1 Score: {f1:.4f}')

# X = pd.concat([X_pos,X_neg], ignore_index=True)

# #### CLEAN DATA ####
# # Add a missing index column to contribute missingness as a feature
# # TODO  it is '.' to indicate missingness
# # TODO Intron_Type change '.' to 'U2'
# # TODO Gene_Regions - change '.' to 'None' (string)

# #neg_train['missing_index'] = neg_train.isnull().any(axis=1).astype(int)
# #pos_train['missing_index'] = pos_train.isnull().any(axis=1).astype(int)

# # Separate features and target variables 
# # X = neg_train.drop(columns=['target'])
# # y = neg_train['target']
# X.Pangolin_Loss = abs(X.Pangolin_Loss)
# X = X.fillna(value=-10)
# print(X)

# #print(X[X.isna().any(axis=1)])

# #X = X[X.columns[~X.columns.isin(['#DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'Pangolin_Gain', 'Pangolin_Loss'])]]

# y = [1] * int(len(X_pos)) + [0] * int(len(X_neg))  
# print(y)

# # Change test size to reflect how big you want you test set to be (0.2-0.4 range)
# # TODO - discuss params, random_state vs stratify
# X_train, X_test, y_train, y_test = train_test_split(
#     X, y, test_size=0.2, random_state=42
# )

# # splice_fit = np.array(X_train[['#DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].max(axis=1)).reshape(-1, 1).tolist() # X_test["SpliceAIMax"].copy()  #X_test.SpliceAIMax.values.reshape(-1, 1)
# # pang_fit = np.array(X_train[['Pangolin_Loss', 'Pangolin_Gain']].max(axis=1)).reshape(-1, 1).tolist() # X_test["PangolinMax"].copy()  #X_test.PangolinMax.values.reshape(-1, 1)

# print("X_train shape : ", X_train.shape) # Train features (without label)
# # print("y_train shape : ", y_train.shape) # Train label of samples
# print("X_test shape : ", X_test.shape) # Test features (without label)
# # print("y_test shape : ", y_test.shape) # Train label of samples
# print(X_train)

# #### TRAIN MODELS ####
# # Initialise classifier 
# # TODO - hyperparameter fine tuning
# rf_classifier = RandomForestClassifier(max_depth=8, random_state=42)
# # Train the model on the training data 
# rf_classifier = rf_classifier.fit(X_train, y_train)
# # Predict labels for the test set features
# y_pred = rf_classifier.predict(X_train)

#### HYPERPARAMETER SELECTION ####
# print("Parameters available : ", rf_classifier.get_params())

# Number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# Number of features to consider at every split
max_features = [None, 'sqrt'] 
# Maximum number of levels in tree
max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]
# Method of selecting samples for training each tree
bootstrap = [True, False]
# Create the random grid
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}
print("Parameter values for testing : ", random_grid)

## Accuracy
print("Accuracy score original: ", accuracy_score(y_test, y_pred))
print("Balanced accuracy score original :" , balanced_accuracy_score(y_test, y_pred))

# exit(0)
#### RECURSIVE FEATURE ELIMINATION ####
# Note: Restarting with a blank model
rfc = regressor
# Remove one feature each step from the model (cross-validated 10 times)
rfc = RFECV(rfc, step=1, cv=10)
rfc = rfc.fit(X_train, y_train)
y_pred_train = rfc.predict(X_train)

print("Feature list : ", rfc.feature_names_in_)
print("Optimal number of features : ", rfc.n_features_)
print("Ranking of features : ", rfc.ranking_)
print("Importance of features ranked #1 : ", rfc.estimator_.feature_importances_)

print(f"{regressor.classes_=}")
exit(0)

# Test the grid of hyperparameters for the best combination
clf_random = RandomizedSearchCV(estimator = rf_classifier, param_distributions = random_grid,
            n_iter = 100, cv = 3, verbose=2, random_state=0, n_jobs = -1, error_score='raise')
clf_grid = clf_random.fit(X_train, y_train)
print("Best parameters : ", clf_grid.best_params_)
# Save best hyperparameters model
best_random_clf = clf_grid.best_estimator_ 
y_pred_train_random = best_random_clf.predict(X_train) 
pickle.dump(best_random_clf, open('RF_train_v1.sav', 'wb'))
print(X_train)
# Precision Recall Curve
PrecisionRecallDisplay.from_estimator(best_random_clf, X_train, y_train, ax = plt.gca(), name="RF CLF")

# Pangolin
pang_X_train = X_train[["Pangolin_Loss", "Pangolin_Gain"]].copy()
p_clf_grid = clf_random.fit(pang_X_train, y_train)
p_best_random_clf = p_clf_grid.best_estimator_ 
pickle.dump(best_random_clf, open('RF_pangolin_train_v1.sav', 'wb'))
PrecisionRecallDisplay.from_estimator(p_best_random_clf, pang_X_train, y_train, ax = plt.gca(), name="RF CLF Pangolin")

# SpliceAI
splice_X_train = X_train[['#DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].copy()
s_clf_grid = clf_random.fit(splice_X_train, y_train)
s_best_random_clf = s_clf_grid.best_estimator_ 
pickle.dump(best_random_clf, open('RF_spliceAI_train_v1.sav', 'wb'))
PrecisionRecallDisplay.from_estimator(s_best_random_clf, splice_X_train, y_train, ax = plt.gca(), name="RF CLF SpliceAI")

# Predict labels for the test set features
# y_pred_test_random = best_random_clf.predict(X_test) 

#### ASSESS MODEL PERFORMANCE ####
## Accuracy
print("Accuracy score original: ", accuracy_score(y_train, y_pred))
print("Balanced accuracy score original :" , balanced_accuracy_score(y_train, y_pred))

print("Accuracy score best hyperparameters: ", accuracy_score(y_train, y_pred_train_random))
print("Balanced accuracy score best hyperparameters:" , balanced_accuracy_score(y_train, y_pred_train_random))

PrecisionRecallDisplay.from_predictions(y_train, splice_fit, ax = plt.gca(), name="X_train.SpliceAIMax")
PrecisionRecallDisplay.from_predictions(y_train, pang_fit, ax = plt.gca(), name="X_train.PangolinMax")
plt.title('Precision-Recall curve')


## Confusion matrix
cm = confusion_matrix(y_train, y_pred_train_random, labels=best_random_clf.classes_)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['Splice Altering', 'Normal'])
disp.plot()
plt.show()

## Classification report
print("Classification report :")
print(classification_report(y_train, y_pred_train_random))

#### RECURSIVE FEATURE ELIMINATION ####
# Note: Restarting with a blank model
rfc = RandomForestClassifier(max_depth=8, random_state=42)
# Remove one feature each step from the model (cross-validated 10 times)
rfc = RFECV(rfc, step=1, cv=10)
rfc = rfc.fit(X_train, y_train)
y_pred_train = rfc.predict(X_train)

print("Feature list : ", rfc.feature_names_in_)
print("Optimal number of features : ", rfc.n_features_)
print("Ranking of features : ", rfc.ranking_)
print("Importance of features ranked #1 : ", rfc.estimator_.feature_importances_)

# %%
