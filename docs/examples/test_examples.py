from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal


def test_l1_regression():
    import l1_regression_gurobipy as gp_impl
    import l1_regression_nupstup as ns_impl

    assert_allclose(gp_impl.y_pred, ns_impl.y_pred)


def test_workforce():
    import workforce_gurobipy as gp_impl
    import workforce_nupstup as ns_impl

    assert_frame_equal(
        gp_impl.assigned_shifts.reset_index(drop=True), ns_impl.assigned_shifts
    )
