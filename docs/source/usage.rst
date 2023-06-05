Usage
=====

The best way to get started with the optimods is to try one out! Check out the
:doc:`gallery` first to find a mod that suits you. Each mod comes complete with
a page of explanatory documentation which provides background on the problem it
is intended to solve. The docs also include runnable example codes and datasets
to get you started.

Note that you may need to install additional dependencies for some examples. The
quick way to ensure you have a Python environment which can execute all the
example snippets in the docs is to run the following::

   python -m pip install gurobi-optimods[examples]

Each mod is designed to be self-contained, with a clean data-in data-out API.
Munge your data into the appropriate format using Python tools you already know,
and run the mod as outlined in its documentation page to get back a solution.

The documentation pages also included details of the mathematical model
underlying the mod, should you be interested to learn more about the
implementation. These details are explained in sections marked "Background". It
is not necessary to read and understand these sections in order to use the mod
effectively. You can also browse :ghsrc:`the mod source <src/gurobi_optimods>`
to find out how the model is implemented in code.

Finally, we welcome contributions of new mods, bug fixes and new features for
existing mods, and improvements to the documentation. This is intended to be a
community project that grows over time to handle a wide range of optimization
use-cases across a different fields. See :doc:`contributing` for more details.
