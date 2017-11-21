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


## Additional Information

The first author -- Dr. Tadamasa Sawada -- presents the highlights of this research in a [YouTube podcast](https://www.youtube.com/watch?v=n5276Nlp16Q).

The abstract of the journal article is reproduced below:

> **Sawada T, Petrov AA.** The divisive normalization model of V1 neurons: a comprehensive comparison of physiological data and model predictions. _J Neurophysiol_ 118: 000–000, 2017. First published August 23, 2017; doi:10.1152/ jn.00821.2016. -- 
> The physiological responses of simple and complex cells in the primary visual cortex (V1) have been studied extensively and modeled at different levels. At the functional level, the divisive normalization model (DNM; Heeger DJ. _Vis Neurosci_ 9: 181–197, 1992) has accounted for a wide range of single-cell recordings in terms of a combination of linear filtering, nonlinear rectification, and divisive normalization. We propose standardizing the formulation of the DNM and implementing it in software that takes static grayscale images as inputs and produces firing rate responses as outputs. We also review a comprehensive suite of 30 empirical phenomena and report a series of simulation experiments that qualitatively replicate dozens of key experiments with a standard parameter set consistent with physiological measurements. This systematic approach identifies novel falsifiable predictions of the DNM. We show how the model simultaneously satisfies the conflicting desiderata of flexibility and falsifiability. Our key idea is that, while adjustable parameters are needed to accommodate the diversity across neurons, they must be fixed for a given individual neuron. This requirement introduces falsifiable constraints when this single neuron is probed with multiple stimuli. We also present mathematical analyses and simulation experiments that explicate some of these constraints.

---
We will try to provide basic support, but cannot make any promises.  
If you discover bugs, please let us know. Thank you!

-- Tadamasa Sawada and Alex Petrov  
2017-11-20
