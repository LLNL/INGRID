INGRID Introduction
===================

INGRID is an interactive grid generator for tokamak edge-plasma modeling capable of handling magnetic-configurations with two x-points in the computational domain.

INGRID analyzes `EFIT <https://gafusion.github.io/OMFIT-source/modules/mod_EFIT.html>`_ data in order to classify the magnetic-topology, and generate the appropriate ``gridue`` formatted file for use in simulation code `UEDGE <https://github.com/LLNL/UEDGE>`_.

Key features of INGRID include:

#. Support for single-null, unbalanced double-null, and snowflake (15, 45, 75, 105, 135, 165) configurations
#. Python based code with GUI and scripting usability
#. UEDGE friendly
#. Portable YAML parameter file driver
#. Modular design pattern for continued development efforts

INGRID was developed at LLNL PLS-FESP by Bryan Garcia (UCSC), Maxim Umansky (LLNL), and Jerome Guterl (GA).