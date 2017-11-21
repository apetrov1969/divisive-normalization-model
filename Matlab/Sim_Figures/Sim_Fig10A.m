%% Simulation Experiments using the standardized Divisive Normalization model
%  Super-saturation effect
%  \beta = 0.0, n_d = 2.35, M = 30;

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt


%% Clean up for debugging
cd(fileparts(mfilename('fullpath'))); % Change directory of this script
clear;
close all;


%% Generate a "model" structure M with default specifications
model_level = 2;
M_SuperSat = DMPL_defaultSpecs([],0,model_level);

M_SuperSat.stim_spec.imageSize_pix = [128,128];
M_SuperSat.EarlyVis_spec.divNorm_exponentDen = 2.35;
M_SuperSat.EarlyVis_spec.divNorm_baselineConst  =0.0;
M_SuperSat.EarlyVis_spec.divNorm_rateScale_sps = 30;
M_SuperSat = DMPL_prepareSpecs(M_SuperSat);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;
neuron_orient_idx = FindClosestIdx(M_SuperSat.EarlyVis_spec.domain_orient_deg,neuron_orient_deg);
neuron_spFreq_idx = FindClosestIdx(M_SuperSat.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);

maxSize_SuperSat = max(abs(M_SuperSat.stim_spec.gridX_deg(:)))*2;


%% Super-saturation effect

    %%  Contrast sensitivity curve
    specs_Contrast.neuron_orient_deg = neuron_orient_deg;
    specs_Contrast.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Contrast.diameter_deg = maxSize_SuperSat;
    data_Contrast = Explore_Contrast(M_SuperSat,specs_Contrast);

    data_Contrast.plot('complex');

    %% Contrast sensitivity by size
    specs_ContBySize.neuron_orient_deg = neuron_orient_deg;
    specs_ContBySize.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_ContBySize = Explore_Contrast_by_Size(M_SuperSat,specs_ContBySize);

    data_ContBySize.plot('complex');
    
    %% Size tuning (contrast)
    specs_SizeByCon.neuron_orient_deg = neuron_orient_deg;
    specs_SizeByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_SizeByCon = Explore_Size_by_Contrast(M_SuperSat,specs_SizeByCon);

    data_SizeByCon.plot('complex');
    
    