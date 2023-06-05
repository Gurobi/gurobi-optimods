License
=======

``gurobi-optimods`` is distributed under the `Apache License 2.0
<https://www.apache.org/licenses/LICENSE-2.0.txt>`_, meaning you can use, read,
and modify the code freely within the license terms. However, it depends on a
commercial package (gurobipy) which requires a license. When installed via pip
or conda, gurobipy ships with a :pypi:`limited trial license <gurobipy>` for the
Gurobi Optimizer which handles models up to a certain size.

The example datasets in the documentation can all be solved using this trial
license, however when using larger input datasets you may encounter the error
"Given data exceeds Gurobi trial license limits". To run the Mod successfully
with your larger dataset, you will need a full license for Gurobi.

- If you are a commercial user: please visit https://www.gurobi.com/free-trial/
  to request a full-featured Evaluation License which will allow you to use
  datasets of unrestricted size in the OptiMods.
- If you are an academic user: Gurobi provides free licenses for faculty
  members, students and staff at recognized academic institutions. In most cases
  you can generate a license yourself through the self-service feature of the
  Gurobi user portal. Please see `How do I obtain a free academic license?
  <https://support.gurobi.com/hc/en-us/articles/360040541251>`_ for further
  instructions.
