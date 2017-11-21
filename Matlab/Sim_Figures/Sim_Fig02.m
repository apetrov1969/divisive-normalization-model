%% Simulation Experiments using the standardized Divisive Normalization model
%  The regular parameter set

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
% Small stimulus size
model_level = 2;
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small = DMPL_prepareSpecs(M_Small);

% Large stimulus size
M_Large = M_Small;
M_Large.stim_spec.imageSize_pix = [128,128];
M_Large = DMPL_prepareSpecs(M_Large);

%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;

maxSize_Small = max(abs(M_Small.stim_spec.gridX_deg(:)))*2;
maxSize_Large = max(abs(M_Large.stim_spec.gridX_deg(:)))*2;


%% Orientation tuning
specs_Orient_CompareModels.neuron_orient_deg = neuron_orient_deg;
specs_Orient_CompareModels.neuron_spFreq_cpd = neuron_spFreq_cpd;
specs_Orient_CompareModels.contrast_pcent = 100;
specs_Orient_CompareModels.diameter_deg = maxSize_Small;
data_Orient_CompareModels = Explore_Orient_CompareModels(M_Small,specs_Orient_CompareModels);

data_Orient_CompareModels.plot();


%% Spatial-frequency tuning
specs_SpFreq_CompareModels.neuron_orient_deg = neuron_orient_deg;
specs_SpFreq_CompareModels.neuron_spFreq_cpd = neuron_spFreq_cpd;
specs_SpFreq_CompareModels.contrast_pcent = 100;
specs_SpFreq_CompareModels.diameter_deg = maxSize_Small;
data_SpFreq_CompareModels = Explore_SpFreq_CompareModels(M_Small,specs_SpFreq_CompareModels);

data_SpFreq_CompareModels.plot();
