import numpy as np
import gurobipy as gp
from gurobipy import GRB

from sklearn import datasets
from sklearn.model_selection import train_test_split

# Load the diabetes dataset
diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

# Split data for fit assessment
X_train, X_test, y_train, y_test = train_test_split(
    diabetes_X, diabetes_y, random_state=42
)

# Metadata about the input data
records, n_features = X_train.shape

# Create model
model = gp.Model()

# Create unbounded variables for each column coefficient, and bound
# magnitudes using additional variables. Keep intercept separate (no
# regularization weight).
intercept = model.addVar(lb=-GRB.INFINITY, name="intercept")
coeffs = np.array(model.addVars(n_features, lb=-GRB.INFINITY, name="coeff").values())
pos_error = np.array(model.addVars(records, name="poserror").values())
neg_error = np.array(model.addVars(records, name="negerror").values())

# Create linear relationship with deviation variables
relation = (X_train * coeffs).sum(axis=1) + intercept + pos_error - neg_error
for i in range(records):
    model.addConstr(relation[i] == y_train[i], name=f"fit_{i}")

# Minimize L1 norm with regularisation term
abs_error = pos_error + neg_error
mean_abs_error = abs_error.sum() / records
model.setObjective(mean_abs_error, sense=GRB.MINIMIZE)

# Optimize to fit coefficients
model.optimize()
intercept = intercept.X
coef = np.array([coeff.X for coeff in coeffs])

# Predict
y_pred = (X_test * coef).sum(axis=1) + intercept
