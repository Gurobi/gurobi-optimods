import unittest

from gurobi_optimods.opf import (
    solve_opf_model,
    generate_opf_solution_figure,
    read_settings_from_file,
    read_coords_from_csv_file,
    read_case_from_file,
    read_case_from_mat_file,
    turn_solution_into_mat_file,
)
from gurobi_optimods.datasets import (
    load_caseopf,
    load_caseopfmat,
    load_caseNYopf,
    load_opfdictcase,
    load_coordsfilepath,
    load_opfgraphicssettings,
    load_opfsettings,
    load_case9solution,
)


class TestOpf(unittest.TestCase):

    numcases = 5
    cases = ["9", "14", "57", "118", "300"]
    # DC test values
    objvals_dc = [5216.026607, 7642.591776, 41006.736942, 125947.881417, 706240.290695]
    Va_dc = [6.177764, 6.283185, 6.171413, 5.817455, -5.520424]
    Pg_dc = [134.377585, 38.032305, 81.931329, 0, 1.2724979]
    Pt_dc = [-56.2622, 69.9608, -18.09715, -102.95381, 24.83]
    # AC test values
    objvals_ac = [5296.686204, 8081.187603]
    Vm_ac = [1.08662, 1.018801]
    Qg_ac = [0.031844, 32.114784]
    Qf_ac = [12.9656, -12.67811]
    # AC relaxation test values
    objvals_acconv = [5296.66532, 8074.9102, 41710.3065, 129338.093, 718613.607]
    Pg_acconv = [89.803524, 194.796114, 142.58252, 24.518669, 0.030902]
    Pt_acconv = [-34.1774, -71.23414, -29.9637, 23.79936, 56.2152]
    # graphics test values
    graphics_9_x = [1129.2, 980.2, 977.6, 1182.8, 480.6, 85.4, 1079.6, 528.0, 0.0]
    graphics_9_y = [
        1066.62,
        132.53,
        220.4,
        0.0,
        777.49,
        569.85,
        1130.71,
        653.08,
        647.86,
    ]

    # Currently, this is just a convenience setting while working on OptiMod
    plot_graphics = False

    # test simple is on purpose the same as test_acopf for now
    # will be removed in final version
    def test_simple(self):
        settings = {
            "doac": True,
            "skipjabr": False,
            "use_ef": True,
            "branchswitching_mip": True,
        }
        # load path to case file
        casefile = load_caseopfmat("9")
        # read case file and return a case dictionary
        case = read_case_from_mat_file(casefile)
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(settings, case, "")
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
        solution, objval = solve_opf_model(settings, case)
        solutionmat, objvalmat = solve_opf_model(settings, casemat)
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertTrue(solutionmat is not None)
        self.assertTrue(objvalmat is not None)

        # objective values should be the same because it's the same data
        self.assertTrue(objval == objvalmat)

        # get path to csv file holding the coordinates for NY
        coordsfile = load_coordsfilepath("nybuses.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        # plot the given solution
        fig = generate_opf_solution_figure({}, case, coords_dict, solution, objval)
        # test a few coordinates
        self.assertLess(abs(fig.data[1].x[0] - 1381.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[0] - 1203.5), 1e-9)
        self.assertLess(abs(fig.data[1].x[-1] - 837.2), 1e-9)
        self.assertLess(abs(fig.data[1].y[-1] - 511.85), 1e-9)
        if self.plot_graphics:
            fig.show()

    # test reading settings and case file from dicts
    def test_opfdicts(self):
        settings = {"branchswitching_mip": True, "doac": True}
        case = load_opfdictcase()
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertLess(abs(solution["f"] - 5296.665647261), 1e-4)
        self.assertLess(abs(solution["bus"][1]["Va"] - 1), 1e-4)
        self.assertLess(abs(solution["gen"][2]["Qg"] - 3.14366), 1e-4)
        self.assertLess(abs(solution["branch"][3]["Pt"] - 56.8647), 1e-4)

    # test reading settings and case file
    def test_settingsfromfile(self):
        settingsfile = load_opfsettings()
        settings = read_settings_from_file(settingsfile)
        case = load_opfdictcase()
        solution, objval = solve_opf_model(settings, case)
        # check whether the solution point looks correct
        self.assertTrue(solution is not None)
        self.assertTrue(objval is not None)
        self.assertLess(abs(solution["f"] - 5296.665647261), 1e-4)
        self.assertLess(abs(solution["bus"][1]["Va"] - 1), 1e-4)
        self.assertLess(abs(solution["gen"][2]["Qg"] - 3.14366), 1e-4)
        self.assertLess(abs(solution["branch"][3]["Pt"] - 56.8647), 1e-4)

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
            self.assertLess(abs(solution["f"] - self.objvals_dc[i]), 1e-4)
            self.assertLess(abs(solutionmat["f"] - self.objvals_dc[i]), 1e-4)
            self.assertLess(abs(solution["bus"][1]["Va"] - self.Va_dc[i]), 1e-4)
            self.assertLess(abs(solutionmat["bus"][1]["Va"] - self.Va_dc[i]), 1e-4)
            self.assertLess(abs(solution["gen"][2]["Pg"] - self.Pg_dc[i]), 1e-4)
            self.assertLess(abs(solutionmat["gen"][2]["Pg"] - self.Pg_dc[i]), 1e-4)
            self.assertLess(abs(solution["branch"][3]["Pt"] - self.Pt_dc[i]), 1e-4)
            self.assertLess(abs(solutionmat["branch"][3]["Pt"] - self.Pt_dc[i]), 1e-4)
            # objective values should be the same because it's the same data and the model is linear
            # without (too much) symmetry
            self.assertLess(abs(objval - objvalmat), 1e-2)

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
            self.assertLess(abs(solution["f"] - self.objvals_ac[i]), 1e-4)
            self.assertLess(abs(solutionmat["f"] - self.objvals_ac[i]), 1e-4)
            self.assertLess(abs(solution["bus"][3]["Vm"] - self.Vm_ac[i]), 1e-4)
            self.assertLess(abs(solutionmat["bus"][3]["Vm"] - self.Vm_ac[i]), 1e-4)
            self.assertLess(abs(solution["gen"][2]["Qg"] - self.Qg_ac[i]), 1e-4)
            self.assertLess(abs(solutionmat["gen"][2]["Qg"] - self.Qg_ac[i]), 1e-4)
            self.assertLess(abs(solution["branch"][1]["Qf"] - self.Qf_ac[i]), 1e-4)
            self.assertLess(abs(solutionmat["branch"][1]["Qf"] - self.Qf_ac[i]), 1e-4)

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
            self.assertLess(abs(objval - self.objvals_acconv[i]), 1)
            self.assertLess(abs(solution["gen"][1]["Pg"] - self.Pg_acconv[i]), 1)
            self.assertLess(abs(solutionmat["gen"][1]["Pg"] - self.Pg_acconv[i]), 1)
            self.assertLess(abs(solution["branch"][2]["Pt"] - self.Pt_acconv[i]), 1)
            self.assertLess(abs(solutionmat["branch"][2]["Pt"] - self.Pt_acconv[i]), 1)

            # objective value should be the same because it's the same data
            # the bigger cases are numerically difficult with large obj vals and we are running Barrier
            # so there can be quite a difference in objective value
            if self.cases[i] == "300":
                # case 300 is quite large and has numerical issues
                self.assertLess(abs(objval - objvalmat), 100)
            else:
                self.assertLess(abs(objval - objvalmat), 10)

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
        self.assertLess(abs(solution["f"] - 5296.716652), 1e-4)
        # TODO finish test for IV
        # print(solution["bus"][3]["Vm"])
        # print(solution["gen"][2]["Qg"])
        # print(solution["branch"][1]["Qf"])

    # test plotting a solution from pre-loaded data
    def test_graphics(self):
        # settings dictionary
        graphics_settings = {}
        # get path to csv file holding the coordinates for case 9
        coordsfile = load_coordsfilepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        # load case dictionary
        case = load_opfdictcase()
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the given solution
        fig = generate_opf_solution_figure(
            graphics_settings, case, coords_dict, solution, objval
        )
        # check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)
        if self.plot_graphics:
            fig.show()

    # test plotting a solution after optimization is performed
    def test_graphics_after_solving(self):
        # load settings dictionary
        optimization_settings = {
            "doac": True,
            "use_ef": True,
            "branchswitching_mip": True,
        }
        # load case dictionary
        case = load_opfdictcase()
        # solve opf model and return a solution and the final objective value
        solution, objval = solve_opf_model(optimization_settings, case)
        # plot the computed solution
        # graphics_settings = {"branchswitching_mip": True}
        graphics_settings = {}
        coordsfile = load_coordsfilepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        fig = generate_opf_solution_figure(
            graphics_settings, case, coords_dict, solution, objval
        )
        # check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)
        if self.plot_graphics:
            fig.show()

    # test plotting a solution while reading graphics settings from a file
    def test_graphics_settings_file(self):
        # load path to settings file
        settingsfile = load_opfgraphicssettings()
        # read settings file and return a settings dictionary
        # set second argument to True because it's a graphics settings file
        graphics_settings = read_settings_from_file(settingsfile, True)
        coordsfile = load_coordsfilepath("case9coords.csv")
        coords_dict = read_coords_from_csv_file(coordsfile)
        # load path to case file
        casefile = load_caseopf("9")
        # read case file and return a case dictionary
        case = read_case_from_file(casefile)
        # load a precomputed solution and objective value
        solution, objval = load_case9solution()
        # plot the computed solution
        fig = generate_opf_solution_figure(
            graphics_settings, case, coords_dict, solution, objval
        )
        # check whether figure coordinates and scaled input coordinates are the same
        for i in range(9):
            self.assertLess(abs(fig.data[1].x[i] - self.graphics_9_x[i]), 1e-9)
            self.assertLess(abs(fig.data[1].y[i] - self.graphics_9_y[i]), 1e-9)
        if self.plot_graphics:
            fig.show()
