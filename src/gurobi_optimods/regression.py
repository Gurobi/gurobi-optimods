import gurobipy as gp
from gurobipy import GRB

from gurobi_optimods.utils import optimod


class LADRegression:
    """Least absolute deviations (L1-norm) regressor"""

    @optimod()
    def fit(self, X_train, y_train, *, create_env):
        """Fit the model to training data.

        :param X_train: Training set feature values
        :type X_train: :class:`np.array`
        :param y_train: Training set output values
        :type y_train: :class:`np.array`
        """

        # Metadata about the input data
        records, self.n_features_in_ = X_train.shape

        # Create model
        with create_env() as env, gp.Model(env=env) as model:

            # Create unbounded variables for each column coefficient, and bound
            # magnitudes using additional variables. Keep intercept separate (no
            # regularization weight).
            intercept = model.addVar(lb=-GRB.INFINITY, name="intercept")
            coeff = model.addMVar(self.n_features_in_, lb=-GRB.INFINITY, name="coeff")
            pos_error = model.addMVar(records, name="pos_error")
            neg_error = model.addMVar(records, name="neg_error")

            # Create linear relationship with deviation variables
            relation = (X_train * coeff).sum(axis=1) + intercept + pos_error - neg_error
            model.addConstr(relation == y_train, name="fit")

            # Minimize L1 norm
            abs_error = pos_error + neg_error
            mean_abs_error = abs_error.sum() / records
            model.setObjective(mean_abs_error, sense=GRB.MINIMIZE)

            # Solve and store results
            model.optimize()
            self.intercept_ = intercept.X
            self.coef_ = coeff.X

    def predict(self, X_test):
        """Predict target value from test data

        :param X_test: Feature data for a new unseen dataset
        :type X_test: :class:`np.array`
        :return: Outputs predicted by the model for the feature data
        :rtype: :class:`np.array`
        """
        return (X_test * self.coef_).sum(axis=1) + self.intercept_
