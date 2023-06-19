Adding a new Mod
================

The goal of the OptiMods is to create and maintain an open-source Python
repository of implemented optimization use cases. Each use case will have clear,
informative, and pretty documentation that explains both how to use it and the
mathematical model behind it for interested readers.

We hope to grow the collection of Mods over time, and welcome your
contributions. This page outlines what makes a good Mod, and what you will need
to consider when developing your contribution.

Proposing a new Mod
-------------------

To propose a new Mod, create an issue in our repository using the 'New Mod
Proposal' template to gather the required details. One of the maintainers will
reach out to you on the issue to discuss the proposed topic and design.

A good Mod:

- solves well-known, well-defined problem from a non-optimization field;
- has a simple interface that shields the user from interacting with gurobipy
  directly;
- is self-contained, and follows a clean data-in data-out style leaning on
  standard packages from the python ecosystem (``numpy``, ``scipy``, and
  ``pandas`` are our first-class citizens);
- is accompanied by clear documentation that provides background information
  for the topic and a formal statement of the problem in the domain language of
  the target user; and
- includes runnable example codes and presentation of results in the
  documentation so users can hit the ground running.

Starting work on a Mod
----------------------

First, carefully read :doc:`contributing` for details on coding standards, and
how to get set up to work on the project. Set up your development environment as
described, and ensure that you can successfully run the tests and build the docs
on your local machine.

Assign yourself (or ask to be assigned) to the relevant Mod proposal issue in
the Github repository. You should then fork `Gurobi/gurobi-optimods
<https://github.com/Gurobi/gurobi-optimods>`_ to your own Github account to
prepare your submission.

Implementation and tests
------------------------

Create ``src/gurobi_optimods/<mod>.py`` where your implementation will live, and
``tests/test_<mod>.py`` where your unit tests will live. A basic implementation
of a Mod takes this form::

    import logging

    import gurobipy as gp

    from gurobi_optimods.utils import optimod

    logger = logging.getLogger(__name__)

    @optimod()
    def my_mod(data, *, create_env):
        """An optimod that solves an important problem

        :param data: Description of argument
        :type data: Type of argument

        ... describe additional arguments ...

        :return: Description of returned result
        :rtype: Type of returned result
        """

        # ... Prepare data ...

        with create_env() as env, gp.Model(env=env) as model:

            # ... Formulate model ...

            model.optimize()

            # ... Extract and post-process the solution ...

            return solution

Mods should be stateless with respect to ``gurobipy`` objects. This means Gurobi
environments and models are created within a Mod function, and closed before the
function returns using context managers. Gurobi environments should be created
by calling ``create_env``. This function is provided to your Mod by the
``@optimod()`` decorator and supplies some necessary parameters to Gurobi to
handle console output and log files consistently across Mods. The standard
parameters ``verbose`` and ``logfile`` are also handled by the decorator and will
be included automatically in the API documentation.

If your Mod needs to produce any output, use the in-built python logging call
``logger.info``.

You should also include your Mod in the :doc:`api` by adding appropriate
`autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_
references to ``docs/source/api.rst``.

Preparing documentation
-----------------------

Create ``docs/source/mods/<mod>.rst`` as the documentation page for your Mod. As
each Mod is different, there is no specific template for this, but please use
the existing Mod pages as a guide.

A reference to your page must also be added to ``docs/source/gallery.rst`` to
include it in the gallery page and toctree when the documentation is built. You
should also add an icon to the gallery card for your Mod.

Your documentation page must contain example codes that new users can
immediately use, with presentation of the results included in the documentation
page.

The implementation of the Mod (description of the mathematical model or
algorithms used) should be hidden from the user at first glance. We use the ``..
dropdown:`` directive for this. Any mathematical and algorithmic details should
be placed in a dropdown and clearly indicate that this is advanced detail the
user does not need to fully understand in order to use the Mod.

Including datasets
------------------

Some of your examples may rely on datasets. These can be packaged with the
optimods to enable users to quickly reproduce the examples in your documentation.

- Any data files should live under a subdirectory
  ``src/gurobi_optimods/data/<mod-name>`` to reduce clutter.
- The ``gurobi_optimods.datasets`` module should implement a function that
  fetches the dataset.

Submitting a pull request
-------------------------

You can submit your pull request at any time in draft mode so the maintainers
are aware your Mod is actively being worked on. This should be a pull request
from a branch on your fork into ``gurobi-optimods/main``. Pull requests include
a checklist of features to ensure your Mod is correctly included in the Python
package and the built documentation.

When your Mod is ready for review, take your PR out of draft mode and request a
review.
