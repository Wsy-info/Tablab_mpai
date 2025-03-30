### Figure 5
### DE

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report

### Heart
### dataset for training the model
adata = anndata.read_h5ad('heart_O_Y_aging.h5ad')
### prepare label
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
y = adata.obs['origin'].values

### split dataset for 70% training set and 30% test set
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify = y, test_size = 0.3, random_state = 2)

### Perform label coding
le = LabelEncoder()
y = le.fit_transform(y)

### Bayesian classifier
nb_model = GaussianNB()
nb_model.fit(X_train, y_train)
nb_predictions = nb_model.predict(X_test)
print("Naive Bayes Classification Report:\n", classification_report(y_test, nb_predictions, target_names=le.classes_))

### dataset for predict
adata_drug = anndata.read_h5ad('heart_aging.h5ad')
### Only retain groups for drug treatments
adata_drug = adata_drug[adata_drug.obs['origin'].isin(['MET', 'NR', 'D+Q', 'SPD', 'MET+NR', 'MET+SPD'])].copy()
x_all = adata_drug.X.toarray() if hasattr(adata_drug.X, "toarray") else adata_drug.X
cell_names = adata_drug.obs.index

### predict the label
y_all = nb_model.predict(x_all)
### print result
results_df = pd.DataFrame({
    'cell_name': cell_names,
    'y_pred': y_all
})
results_df.to_csv('prediction_heart_aging_bys.csv', index = False)