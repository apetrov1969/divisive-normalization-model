%% Simulation Experiments using the standardized Divisive Normalization model

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

%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;

maxSize_Small = max(abs(M_Small.stim_spec.gridX_deg(:)))*2;

    
%% Contrast Response of Plaid

%% Standard parameter set
specs_PlaidCont.neuron_orient_deg = neuron_orient_deg;
specs_PlaidCont.neuron_spFreq_cpd = neuron_spFreq_cpd;
specs_PlaidCont.diameter_deg = maxSize_Small;
data_PlaidCont = Explore_PlaidContrast(M_Small,specs_PlaidCont);

data_PlaidCont.plot('complex');


%% Modified parameter set A (h_\Theta = 90)
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small.EarlyVis_spec.fwmh_pool_surr_orient_deg = 90; %60
M_Small = DMPL_prepareSpecs(M_Small);

data_PlaidCont = Explore_PlaidContrast(M_Small,specs_PlaidCont);
data_PlaidCont.plot('complex');


%% Modified parameter set B (n_n = n_d = 10)
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small.EarlyVis_spec.divNorm_exponentNum = 10; % 2.0
M_Small.EarlyVis_spec.divNorm_exponentDen = 10; % 2.0
M_Small = DMPL_prepareSpecs(M_Small);

data_PlaidCont = Explore_PlaidContrast(M_Small,specs_PlaidCont);
data_PlaidCont.plot('complex');


%% Modified parameter set C (h_\Theta = 90, n_n = n_d = 10)
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small.EarlyVis_spec.fwmh_pool_surr_orient_deg = 90; %60
M_Small.EarlyVis_spec.divNorm_exponentNum = 10; % 2.0
M_Small.EarlyVis_spec.divNorm_exponentDen = 10; % 2.0
M_Small = DMPL_prepareSpecs(M_Small);

data_PlaidCont = Explore_PlaidContrast(M_Small,specs_PlaidCont);
data_PlaidCont.plot('complex');



