Welcome to Gurobi OptiMods's documentation!
===========================================

**Gurobi OptiMods**: nice APIs for common optimization tasks.

``gurobi-optimods`` is an open-source Python repository of implemented
optimization use cases, each with clear, informative, and pretty documentation
that explains how to use it and the mathematical model behind it.

Getting Started
---------------

The package is a collection of independent 'mods'. Each mod is intended to be
immediately applicable to real use-cases. However, we expect that for many
practical applications users will need to understand and extend the
implementation of a mod to tailor it to their use-case. Read the :doc:`usage`
section first for an overview of the design and use-case for the OptiMods.

Check out the :doc:`mods/index` for a quick overview of the current set of
implemented mods. We welcome contributions of use-cases you are interested in,
or extensions to existing mods. See the :doc:`contributing` for more
information on how to get involved in the project.

.. note::

   This project is under active development, Gurobi-OptiMods has its
   documentation hosted on Read the Docs.

Contents
--------

.. toctree::
   :maxdepth: 1
   :caption: Start

   installation
   usage

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   mods/index
   extending
   contributing

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api
   license
   contact
