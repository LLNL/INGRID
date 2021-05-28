# INGRID

Tokamak edge plasma grid generator for magnetic-topologies with up to two x-points anywhere in the domain.

## Description

INGRID (Interactive Grid) is a Python based tokamak edge plasma grid generator capable of automatic generation of grids for magnetic-topologies with up to two x-points anywhere in the domain. The code generates block-structured grids in three steps:

    - Analysis of MHD reconscruction code such as EFIT to identify the embedded topology.
    - Generation of the appropriate Patch map which appropriately models the domain (blocks).
    - Generation of the grid ready for export.
    
INGRID can be utilized in a GUI mode and noninteractively via Python scripts. Both modes of operation support the use of the YAML formatted parameter file.

All documentation and tutorials pertaining to INGRID are hosted on Read The Docs, and can be found [here](https://ingrid.readthedocs.io/en/latest/).

## Obtaining INGRID

Installation instructions can be found [here](https://ingrid.readthedocs.io/en/latest/installation.html) on Read The Docs.

## Getting Started
Instructions for starting an INGRID session and generating grids are [here](https://ingrid.readthedocs.io/en/latest/getting_started.html) on Read The Docs.

## Development Information

INGRID was developed at Lawrence Livermore National Laboratory by Bryan Garcia (UCSC), Jerome Guterl (GA), Joey Watkins (BYU), and Maxim Umansky (LLNL). Maintained by Bryan Garcia and Maxim Umansky.