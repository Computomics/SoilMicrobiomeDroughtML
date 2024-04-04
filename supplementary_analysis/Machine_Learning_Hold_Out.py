import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
import pickle

"""
This script performs nested cross-validation with hyperparameter tuning 
on the grass-drought dataset after removing 20% of the samples randomly 
while still preserving the ratio between 'Drought' and 'Control' samples 
and using the samples that were removed as an independent hold-out dataset. 
Both datasets are preprocessed exactly the same as the initial grass-drought 
dataset and feature tables were created where the reduced dataset contained 
448 samples and the hold-out dataset contained 112 samples. 
Performance was tested for the nested cross-validation and the independent 
hold-out dataset.
"""

ranks = ['phylum', 'class', 'order', 'family', 'genus']
ml_model = 'RFC'
    
for rank in ranks:

    df = pd.read_csv(f'feature_tbl_{rank}.csv')
    df = df.drop(df.columns[0], axis=1)

    X, y = df.drop('Watering_Regm', axis=1), df['Watering_Regm']

    cv_outer = KFold(n_splits=5, shuffle=True, random_state=42)
    

    # performance reports
    accuracy_results = list()
    f1_results = list()
    precision_results = list()
    recall_results = list()
    auc_results = list()

    df = df.drop('Watering_Regm', axis = 1)
    feature_imp_coeff = pd.DataFrame(X.columns, columns=['features'])
    enriched_all = pd.DataFrame(X.columns, columns = ['Taxa'])
    list_shap_values = list()
    list_test_sets = list()
    idx = 0

    best_params_per_fold = {}

    for train_ix, test_ix in cv_outer.split(X):
        
        # split data
        X_train, X_test = X.iloc[train_ix], X.iloc[test_ix]
        y_train, y_test = y.iloc[train_ix], y.iloc[test_ix]

        # configure the cross-validation procedure
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)

        # define the model
        model = RandomForestClassifier(random_state=42)

        # define search space
        space = dict()
        space['n_estimators'] = [100, 200, 500, 750, 1000, 1500, 2000]
        space['max_depth'] = [int(x) for x in np.linspace(10, 110, num = 11)]
        space['min_samples_leaf'] = [1, 2, 4]
        space['min_samples_split'] = [2, 5, 10]
        space['max_features'] = ['sqrt', 'log2']
            
        # define search
        search = GridSearchCV(model, space, scoring='accuracy', n_jobs=80, cv=cv_inner, refit=True)
        # execute search
        result = search.fit(X_train, y_train)

        # get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        
        # get params of best model
        best_params_for_fold = result.best_params_

        # evaluate model on the hold out dataset
        y_pred = best_model.predict(X_test)
        # evaluate the model
        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred)
        rec = recall_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_pred)

        # store the result
        accuracy_results.append(acc)
        f1_results.append(f1)
        precision_results.append(prec)
        recall_results.append(rec)
        auc_results.append(auc)

        # save best estimator
        filename = f'{rank}_{ml_model}_{idx+1}.sav'
        pickle.dump(best_model, open(filename, 'wb'))
        
        idx += 1

    # summarize the estimated performance of the model
    with open(f'performance_{rank}_{ml_model}.txt', 'w') as file:

        file.write('Accuracy: %.3f (%.3f), F1: %.3f (%.3f), Precision: %.3f (%.3f), Recall: %.3f (%.3f), AUC: %.3f (%.3f)\n' % (
            np.mean(accuracy_results), np.std(accuracy_results),
            np.mean(f1_results), np.std(f1_results),
            np.mean(precision_results), np.std(precision_results),
            np.mean(recall_results), np.std(recall_results),
            np.mean(auc_results), np.std(auc_results)
            
        ))


# Predicting on the unseen hold-out dataset of the grass-drought dataset

for rank in ranks:

    # Reading in test_dataset
    df_test = pd.read_csv(f'hold_out_feature_tbl_{rank}.csv')
    df_test = df_test.drop(df_test.columns[0], axis=1)
    X_test, y_test = df_test.drop('Watering_Regm', axis=1), df_test['Watering_Regm']

    # performance reports
    accuracy_results = list()
    f1_results = list()
    precision_results = list()
    recall_results = list()
    auc_results = list()
    
    for i in range(5):
        
        model = pickle.load(open(f'{rank}_{ml_model}_{i+1}.sav', 'rb'))
        y_pred = model.predict(X_test)
        
        accuracy_results.append(accuracy_score(y_test, y_pred))
        f1_results.append(f1_score(y_test, y_pred))
        precision_results.append(precision_score(y_test, y_pred))
        recall_results.append(recall_score(y_test, y_pred))
        
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        auc_results.append(roc_auc_score(y_test, y_pred_proba))
        
    with open(f'hold_out_performance_{rank}_{ml_model}.txt', 'w') as file:
        
        file.write('Accuracy: %.3f (%.3f), F1: %.3f (%.3f), Precision: %.3f (%.3f), Recall: %.3f (%.3f), AUC: %.3f (%.3f)\n' % (
            np.mean(accuracy_results), np.std(accuracy_results),
            np.mean(f1_results), np.std(f1_results),
            np.mean(precision_results), np.std(precision_results),
            np.mean(recall_results), np.std(recall_results),
            np.mean(auc_results), np.std(auc_results)
            
    ))