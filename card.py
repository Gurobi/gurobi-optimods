### First snippet
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from gurobi_optimods.regression import CardinalityConstrainedRegression, LADRegression

# Load the diabetes dataset
diabetes = datasets.load_diabetes()

# Split data for fit assessment
X_train, X_test, y_train, y_test = train_test_split(
    diabetes["data"], diabetes["target"], random_state=982, test_size=0.1
)

# Fit model and obtain predictions
ccr = CardinalityConstrainedRegression(k=10)
ccr.fit(X_train, y_train, silent=False)
y_pred = ccr.predict(X_test)

### Second snippet
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso

lasso = Lasso(alpha=0.1)
lasso.fit(X_train, y_train)
coefficients = pd.DataFrame(
    data={"Lasso": lasso.coef_, "CCR": ccr.coef_}, index=diabetes["feature_names"]
)

fig = plt.figure(figsize=(8, 4))
coefficients.plot.bar(ax=plt.gca())
fig.tight_layout()

### Hyperparameter optimization for CCR

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso, Ridge, LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error

models = [
    {"k": k, "regressor": CardinalityConstrainedRegression(k=k)} for k in range(2, 10)
]

for model in models:
    model["regressor"].fit(X_train, y_train, silent=True)
    y_pred = model["regressor"].predict(X_train)
    model["mae-train"] = mean_absolute_error(y_pred, y_train)
    model["mse-train"] = mean_squared_error(y_pred, y_train)
    y_pred = model["regressor"].predict(X_test)
    model["mae-test"] = mean_absolute_error(y_pred, y_test)
    model["mse-test"] = mean_squared_error(y_pred, y_test)

pd.DataFrame(models).plot.line(x="k", y=["mae-test", "mae-train"])

### Hyperparameter optimization for Lasso

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error

models = [
    {"alpha": alpha, "regressor": Lasso(alpha=alpha)}
    for alpha in np.linspace(0.001, 1.0)
]

for model in models:
    model["regressor"].fit(X_train, y_train)
    y_pred = model["regressor"].predict(X_train)
    model["mae-train"] = mean_absolute_error(y_pred, y_train)
    model["mse-train"] = mean_squared_error(y_pred, y_train)
    y_pred = model["regressor"].predict(X_test)
    model["mae-test"] = mean_absolute_error(y_pred, y_test)
    model["mse-test"] = mean_squared_error(y_pred, y_test)

pd.DataFrame(models).plot.line(x="alpha", y=["mae-test", "mae-train"])
