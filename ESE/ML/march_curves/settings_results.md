# SpliceAI (0.97)
```python
regressor = RandomForestClassifier(n_estimators=30,
                                   random_state=42,
                                   max_depth=5,
                                   min_samples_split=5,
                                   min_samples_leaf=5,
                                   criterion='entropy',
                                   class_weight='balanced_subsample', n_jobs=-1)
```

```
Training:
F1 Score: 0.928
Accuracy score (training): 0.924
Balanced accuracy score (training): 0.927
Validation:
F1 Score: 0.926
Accuracy score (test): 0.923
Balanced accuracy score (test): 0.925
```

# Pangolin (0.97)
```python
regressor = RandomForestClassifier(n_estimators=30,
                                   random_state=42,
                                   max_depth=5,
                                   min_samples_split=5,
                                   min_samples_leaf=5,
                                   criterion='entropy',
                                   class_weight='balanced_subsample', n_jobs=-1)
```

```
Training:
F1 Score: 0.937
Accuracy score (training): 0.933
Balanced accuracy score (training): 0.936
Validation:
F1 Score: 0.932
Accuracy score (test): 0.928
Balanced accuracy score (test): 0.930
```
---
```python
regressor = LogisticRegression(class_weight='balanced')
```

```
Training:
F1 Score: 0.918
Accuracy score (training): 0.916
Balanced accuracy score (training): 0.922
Validation:
F1 Score: 0.916
Accuracy score (test): 0.913
Balanced accuracy score (test): 0.918
```
---
```python
regressor = SVC(probability=True, class_weight='balanced',
                **{'kernel': 'rbf', 'gamma': np.float64(1.0), 'degree': 2, 'C': np.float64(35.93813663804626)})
    # grid search
    param_dist = {
        'C': np.logspace(-2, 2, 10),  # e.g., 0.01 to 100
        'kernel': ['linear', 'rbf', 'poly'],
        'gamma': ['scale', 'auto'] + list(np.logspace(-3, 0, 5)),
        'degree': [2]  # relevant only if kernel='poly'
    }

    random_search = RandomizedSearchCV(
        estimator=regressor,
        param_distributions=param_dist,
        n_iter=20,   # number of random combinations to try
        scoring='f1',
        cv=5,
        random_state=42,
        n_jobs=-1,
        verbose=2
    )
```

```
Training:
F1 Score: 0.933
Accuracy score (training): 0.930
Balanced accuracy score (training): 0.933
Validation:
F1 Score: 0.932
Accuracy score (test): 0.928
Balanced accuracy score (test): 0.930
```

# ESE (0.72)

```python
regressor = RandomForestClassifier(n_estimators=30,
                                   random_state=42,
                                   max_depth=5,
                                   min_samples_split=5,
                                   min_samples_leaf=5,
                                   criterion='entropy',
                                   class_weight='balanced_subsample', n_jobs=-1)
```

```
Training:
F1 Score: 0.595
Accuracy score (training): 0.641
Balanced accuracy score (training): 0.658
Validation:
F1 Score: 0.601
Accuracy score (test): 0.648
Balanced accuracy score (test): 0.661
```
---
```python
# Grid search
regressor = SVC(probability=True, class_weight='balanced',
                **{'kernel': 'rbf', 'gamma': np.float64(1.0), 'C': np.float64(0.01)})
```

```
Training:
F1 Score: 0.700
Accuracy score (training): 0.576
Balanced accuracy score (training): 0.541
Validation:
F1 Score: 0.693
Accuracy score (test): 0.565
Balanced accuracy score (test): 0.535
```
---
```python
regressor = LogisticRegression(class_weight='balanced')
```

```
Training:
F1 Score: 0.627
Accuracy score (training): 0.641
Balanced accuracy score (training): 0.651
Validation:
F1 Score: 0.625
Accuracy score (test): 0.638
Balanced accuracy score (test): 0.645
```

# mmSplice (0.94)

Random Forest 
```
Training:
F1 Score: 0.878
Accuracy score (training): 0.878
Balanced accuracy score (training): 0.884
Validation:
F1 Score: 0.872
Accuracy score (test): 0.871
Balanced accuracy score (test): 0.879
```

# spliceogen (0.95)

dots to zeros

Random Forest
```
Training:
F1 Score: 0.895
Accuracy score (training): 0.891
Balanced accuracy score (training): 0.895
Validation:
F1 Score: 0.900
Accuracy score (test): 0.896
Balanced accuracy score (test): 0.899
```

# spip (0.94)

Random Forest
```
Training:
F1 Score: 0.880
Accuracy score (training): 0.876
Balanced accuracy score (training): 0.880
Validation:
F1 Score: 0.879
Accuracy score (test): 0.874
Balanced accuracy score (test): 0.878
```
