from sklearn import datasets
from sklearn.model_selection import train_test_split

from nupstup.regression import GurobiL1Regression

# Load the diabetes dataset
diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

# Split data for fit assessment
X_train, X_test, y_train, y_test = train_test_split(
    diabetes_X, diabetes_y, random_state=42
)

# Create and fit parameterised model
reg = GurobiL1Regression(alpha=0.01)
reg.fit(X_train, y_train)
y_pred = reg.predict(X_test)
