import numpy as np
import gurobipy as gp
from gurobipy import GRB


class GurobiL1Regression:
    """ L1-norm regressor with Lasso regularization """

    def __init__(self, alpha):
        self.alpha = alpha

    def fit(self, X_train, y_train):
        """Fit the model to training data.
        """

        # Metadata about the input data
        records, self.n_features_in_ = X_train.shape

        # Create model
        model = gp.Model()

        # Create unbounded variables for each column coefficient, and bound
        # magnitudes using additional variables. Keep intercept separate (no
        # regularization weight).
        intercept = model.addVar(lb=-GRB.INFINITY, name="intercept")
        coeffs = np.array(
            model.addVars(self.n_features_in_, lb=-GRB.INFINITY, name="coeff").values()
        )
        abscoeff = np.array(
            model.addVars(self.n_features_in_, name="abscoeff").values()
        )
        pos_error = np.array(model.addVars(records, name="poserror").values())
        neg_error = np.array(model.addVars(records, name="negerror").values())

        for i in range(self.n_features_in_):
            model.addConstr(coeffs[i] <= abscoeff[i], name=f"poscoeff_{i}")
            model.addConstr(coeffs[i] >= -abscoeff[i], name=f"negcoeff_{i}")

        # Create linear relationship with deviation variables
        relation = (X_train * coeffs).sum(axis=1) + intercept + pos_error - neg_error
        for i in range(records):
            model.addConstr(relation[i] == y_train[i], name=f"fit_{i}")

        # Minimize L1 norm with regularisation term
        abs_error = pos_error + neg_error
        mean_abs_error = abs_error.sum() / records
        regularization = abscoeff.sum()
        model.setObjective(
            mean_abs_error + self.alpha * regularization, sense=GRB.MINIMIZE
        )

        # Optimize
        model.optimize()

        self.intercept_ = intercept.X
        self.coef_ = np.array([coeff.X for coeff in coeffs])

    def predict(self, X_test):
        """Predict target value from test data
        """
        return (X_test * self.coef_).sum(axis=1) + self.intercept_
