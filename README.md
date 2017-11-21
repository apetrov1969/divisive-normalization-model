# Divisive Normalization Model

This repository contains Matlab software associated with the following publication:

**Sawada, T. & Petrov, A. A. (2017). The divisive normalization model of V1 neurons: A comprehensive comparison of physiological data and model predictions. _Journal of Neurophysiology, 118_, xxx-xxx. [doi:10.1152/jn.00821.2016](https://doi.org/10.1152/jn.00821.2016).**

Copyright (C) 2013-2017 Laboratory for Cognitive Modeling and Computational Cognitive Neuroscience at the Ohio State University, <http://cogmod.osu.edu>.


The software in this repository is released under [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html). See the `LICENSE` file for details.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


## Matlab Files

The `Matlab` directory contains the file `Installer_Matlab.m` that adds the relevant files to the Matlab path.  There are the following subdirectories:

* `StandardDivisiveNormalization` -- Matlab functions that implement the Standard Divisive Normalization Model (DNM).  Note that the `DMPL_` prefix is synonymous with DNM.
* `Utility` -- Short utility functions.
* `Sim_EachExperiment` -- Scripts that explore various phenomena (cf. Table 1 in the paper). Only some of these simulation results are described in the journal article.
* `Sim_Figures` -- Scripts that generate the simulation-related figures of the journal article.
* `SampleImages` -- Bitmap (`.bmp`) images of various gratings that can be used as inputs to the model.


## Getting Started

The `UserManual` directory contains an all-too-short document `UserManual.docx` that explains how to run the model.  The associated Matlab scripts illustrate the basic functionality.

We would try to provide basic support, but cannot make any promises.  
If you discover any bugs, please let us know. Thank you!

-- Tadamasa Sawada and Alex Petrov

November 2017
