import unittest

from gurobi_optimods.opf import (
    solve_opf_model,
    plot_opf_solution,
    read_settings_from_file,
    read_case_from_file,
    read_case_from_mat_file,
)
from gurobi_optimods.datasets import (
    load_caseopf,
    load_caseopfmat,
    load_caseNYopf,
    load_opfdictcase,
    load_opfgraphicssettings,
    load_opfsettings,
    load_case9solution,
)


class TestOpf(unittest.TestCase):

    numcases = 5
    cases = ["9", "14", "57", "118", "300"]
    # DC test values
    numvars_dc = [50, 95, 363, 850, 1904]
    objvals_dc = [5216.026607, 7642.591776, 41006.736942, 125947.881417, 706240.290695]
    lincost_dc = [688.133508, 5179.999999, 29931.879752, 84840, 470847.849469]
    theta_9_dc = [6.083, 6.008502, 6.177336, 6.12929, -5.579928]
    # AC test values
    numvars_ac = [107, 202]
    objvals_ac = [5296.686204, 8081.187603]
    lincost_ac = [704.365640, 6108.954838]
    e_1_ac = [1.0999999, -1.0599999]

    # AC relaxation test values
    numvars_acconv = [89, 174, 667, 1580, 3506]
    objvals_acconv = [5296.66532, 8074.9102, 41710.3065, 129338.093, 718613.607]
    lincost_acconv = [704.38492, 6095.638486, 30498.743097, 93490.138573, 485687.552521]
    IQ_1_acconv = [0.129955, -0.0564129, 0.319624, -0.1200001, -0.49]

    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        settings = {"doac": True, "skipjabr": False, "use_ef": True}
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case, "OPF.log")
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)

    # test a real data set for New York
    def test_NY(self):
        settings = {"dodc": True}
        # load path to case file
        casefile, casemat = load_caseNYopf()
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        casemat = read_case_from_mat_file(casemat)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case, "OPF.log")
        solutionmat, objvalmat = solve_opf_model(settings, casemat, "OPF.log")
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(solutionmat is not None)
        self.assertTrue(objvalmat is not None)

        # solutions should be the same because it's the same data
        self.assertTrue(objval == objvalmat)
        for s in solution:
            self.assertTrue(solution[s] == solutionmat[s])

    # test reading settings and case file from dicts
    def test_opfdicts(self):
        settings = {"branchswitching_mip": True, "doac": True}
        case = load_opfdictcase()
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.399018), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-4)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.127732), 1e-4)

    # test reading settings and case file
    def test_settingsfromfile(self):
        settingsfile = load_opfsettings()
        settings = read_settings_from_file(settingsfile)
        case = load_opfdictcase()
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 134)
        self.assertLess(abs(objval - 5296.665647261), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 704.399018), 1e-4)
        self.assertTrue("twinP_1_1_4" in solution.keys())
        self.assertLess(abs(solution["twinP_1_1_4"]), 1e-4)
        self.assertTrue("Q_9_4_9" in solution.keys())
        self.assertLess(abs(solution["Q_9_4_9"] - 0.127732), 1e-4)

    # test DC formulation
    def test_dcopf(self):
        # set settings
        settings = {"branchswitching_mip": True, "dodc": True}

        for i in range(self.numcases):
            # load path to case file in .m and .mat format
            casefile = load_caseopf(self.cases[i])
            casefilemat = load_caseopfmat(self.cases[i])
            # read case file in .m and .mat format and return a case dictionary
            case = read_case_from_file(casefile)
            casemat = read_case_from_mat_file(casefilemat)
            # solve opf models and return a solution and the final objective value
            solution, objval = solve_opf_model(settings, case)
            solutionmat, objvalmat = solve_opf_model(settings, casemat)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(objval is not None)
            self.assertTrue(solutionmat is not None)
            self.assertTrue(objvalmat is not None)
            self.assertTrue(len(solution) == self.numvars_dc[i])
            self.assertLess(abs(objval - self.objvals_dc[i]), 1e-4)
            self.assertTrue("lincost" in solution.keys())
            self.assertLess(abs(solution["lincost"] - self.lincost_dc[i]), 1e-4)
            self.assertTrue("theta_9" in solution.keys())
            self.assertLess(abs(solution["theta_9"] - self.theta_9_dc[i]), 1e-4)

            # solutions should be the same because it's the same data and the model is linear
            # without (too much) symmetry
            self.assertLess(abs(objval - objvalmat), 1e-2)
            for s in solution:
                self.assertLess(abs(solution[s] - solutionmat[s]), 1e-2)

    # test AC formulation
    def test_acopf(self):
        settings = {"doac": True, "use_ef": True}

        for i in range(2):
            # load path to case file in .m and .mat format
            casefile = load_caseopf(self.cases[i])
            casefilemat = load_caseopfmat(self.cases[i])
            # read case file in .m and .mat format and return a case dictionary
            case = read_case_from_file(casefile)
            casemat = read_case_from_mat_file(casefilemat)
            # solve opf models and return a solution and the final objective value
            solution, objval = solve_opf_model(settings, case)
            solutionmat, objvalmat = solve_opf_model(settings, casemat)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(objval is not None)
            self.assertTrue(solutionmat is not None)
            self.assertTrue(objvalmat is not None)
            self.assertTrue(len(solution) == self.numvars_ac[i])
            self.assertLess(abs(objval - self.objvals_ac[i]), 1)
            self.assertTrue("lincost" in solution.keys())
            self.assertLess(abs(solution["lincost"] - self.lincost_ac[i]), 1e-4)
            self.assertTrue("e_1" in solution.keys())
            self.assertLess(abs(solution["e_1"] - self.e_1_ac[i]), 1e-4)

            # objective value should be the same because it's the same data
            self.assertLess(abs(objval - objvalmat), 1e-2)
            # cannot check solution point because it can be symmetric

    # test AC formulation relaxation
    def test_acopfconvex(self):
        settings = {"doac": True}

        for i in range(self.numcases):
            # load path to case file in .m and .mat format
            casefile = load_caseopf(self.cases[i])
            casefilemat = load_caseopfmat(self.cases[i])
            # read case file in .m and .mat format and return a case dictionary
            case = read_case_from_file(casefile)
            casemat = read_case_from_mat_file(casefilemat)
            # solve opf models and return a solution and the final objective value
            solution, objval = solve_opf_model(settings, case)
            solutionmat, objvalmat = solve_opf_model(settings, casemat)
            # check whether the solution point looks correct
            self.assertTrue(solution is not None)
            self.assertTrue(objval is not None)
            self.assertTrue(solutionmat is not None)
            self.assertTrue(objvalmat is not None)
            self.assertTrue(len(solution) == self.numvars_acconv[i])
            self.assertLess(abs(objval - self.objvals_acconv[i]), 1)
            self.assertTrue("lincost" in solution.keys())
            self.assertLess(abs(solution["lincost"] - self.lincost_acconv[i]), 1e-4)
            self.assertTrue("IQ_1" in solution.keys())
            self.assertLess(abs(solution["IQ_1"] - self.IQ_1_acconv[i]), 1e-4)

            # objective value should be the same because it's the same data
            # the bigger cases are numerically difficult with large obj vals and we are running Barrier
            # so there can be quite a difference in objective value
            if self.cases[i] == "300":
                # case 300 is quite large and has numerical issues
                self.assertLess(abs(objval - objvalmat), 100)
            else:
                self.assertLess(abs(objval - objvalmat), 1)

    # test IV formulation
    def test_ivopf(self):
        # set settings
        settings = {"doiv": True, "ivtype": "aggressive"}
        # load path to case file
        # currently all other cases take very long in IV formulation
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution points looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(len(solution) == 62)
        self.assertLess(abs(objval - 5296.716652), 1e-4)
        self.assertTrue("lincost" in solution.keys())
        self.assertLess(abs(solution["lincost"] - 703.019340), 1e-4)
        self.assertTrue("f_9" in solution.keys())
        self.assertLess(abs(solution["f_9"] + 0.085954), 1e-4)
        self.assertTrue("P_9_4_9" in solution.keys())
        self.assertLess(abs(solution["P_9_4_9"] - 0.543771), 1e-4)

    def test_graphics(self):
        # settings dictionary
        settings_graphics_dict = {"branchswitching_mip": True}
        # load case dictionary
        case = load_opfdictcase()
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the given solution
        plot_opf_solution(settings_graphics_dict, case, solution, objval)

    def test_graphics_after_solving(self):
        # load settings dictionary
        settings_optimization_dict = {"branchswitching_mip": True, "doac": True}
        # load case dictionary
        case = load_opfdictcase()
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings_optimization_dict, case)
        # plot the computed solution
        settings_graphics_dict = {"branchswitching_mip": True}
        plot_opf_solution(settings_graphics_dict, case, solution, objval)

    def test_graphics_settings_file(self):
        # load path to settings file
        settingsfile = load_opfgraphicssettings()
        # read settings file and return a settings dictionary
        settings_graphics_dict = read_settings_from_file(settingsfile, True)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the computed solution
        plot_opf_solution(settings_graphics_dict, case, solution, objval)
