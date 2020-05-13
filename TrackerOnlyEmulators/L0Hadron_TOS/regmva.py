
import random

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import classification_report, roc_auc_score



from root_numpy import root2array, rec2array


branch_names = """Dplus_PT,Dplus_P,Kminus_L0Calo_HCAL_realET,piplus_L0Calo_HCAL_realET,piplus0_L0Calo_HCAL_realET,pipi_diff,piK_diff,pi0K_diff""".split(",")
_names = [c.strip() for c in branch_names]
branch_names = (b.replace(" ", "_") for b in branch_names)
branch_names = list(b.replace("-", "_") for b in branch_names)
variables = root2array("regmva.root",
                    "DecayTree", 
                    branch_names)
X = rec2array(variables)

real_triggerET = root2array("regmva.root",
                    "DecayTree", 
                    "DTrigger_real")

emulated_triggerET = root2array("regmva.root",
                    "DecayTree", 
                    "DTrigger_emulated")

# for sklearn data is usually organised
# into one 2D array of shape (n_samples x n_features)
# containing all the data and one array of categories
# of length n_samples



y = real_triggerET-emulated_triggerET


import pandas.core.common as com
from pandas.core.index import Index

from pandas import plotting
from pandas.plotting import scatter_matrix


from sklearn.model_selection import train_test_split

X_dev,X_eval, y_dev,y_eval = train_test_split(X, y,
                                              test_size=0.33, random_state=42)
X_train,X_test, y_train,y_test = train_test_split(X_dev, y_dev,
                                                  test_size=0.33, random_state=492)


from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
rng = np.random.RandomState(1)
bdt = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
                          n_estimators=300, random_state=rng)

bdt.fit(X_train, y_train)
print('BDT fitted')

import pickle
pickle.dump(bdt,open("bdt.joblib","wb"))

bdt = pickle.load(open("bdt.joblib","rb"))


y_predicted = bdt.predict(X)
y_predicted.dtype = [('Trigger_correction', 'float64')]

print(type(y_predicted))
print(len(y_predicted))
from root_numpy import array2root



array2root(y_predicted, "RegBDT_predicted.root", "BDToutput", mode='recreate')




