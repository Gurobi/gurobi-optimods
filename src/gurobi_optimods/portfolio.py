import gurobipy as gp

# In general the public API should be a single class or function. Go with
# whatever makes the most sense for this mod.

class Portfolio:
    def __init__(
        self,
        H,
        Sigma,
        mu,
        current_alloc=None,
        maxnum_transactions=20,
        maxnum_allocations=50,
        min_invest_long=0.05,
        min_invest_short=0.01,
        max_total_short=0.1,
    ) -> None:

        self.factor_matrix = H  # H.T@H is the variance-covariance matrix
        self.covariance = Sigma  # Sigma is the variance-covariance matrix
        self.mu = mu  # estimated first moments of return function
        self.current_alloc = current_alloc  # existing allocation
        self.maxnum_transactions = maxnum_transactions  # No more than 20 trades
        self.maxnum_allocations = maxnum_allocations  # No more than 50 allocations at a time (e.g., online problem)
        self.min_invest_long = min_invest_long  # For long allocations, need at least 5% of total investment
        self.min_invest_short = min_invest_short  # For short allocations, need at least 1% of total investment
        self.max_total_short = max_total_short  # Maximum 10% short selling


def solve_mod(data):
    pass


def minimize_risk(data):
    """
    A sphinx-compatible docstring

    :param data1: Data structure for first argument
    :type data1: pd.DataFrame
    """

    # min x.T @ H.T @ H @ x - mu @ x
    # s.t. something

    with gp.Env() as env, gp.Model(env=env) as model:
        # build model
        model.optimize()
        # post-process and return solution
        return


def maximize_return(data):
    pass
    # max  mu @ x
    # s.t. something
    #     x.T @ H.T @ H @ x <= omega^2
