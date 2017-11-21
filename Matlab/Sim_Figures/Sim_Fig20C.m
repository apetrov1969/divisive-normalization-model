%% Simulation Experiments using the standardized Divisive Normalization model
%  Effect of the surround suppression on contrast sensitivity curve
%  n_n = 3.0, n_d = 2.8

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
M_nn30nd28 = DMPL_defaultSpecs([],0,model_level);
M_nn30nd28.stim_spec.imageSize_pix = [128,128];
M_nn30nd28.EarlyVis_spec.divNorm_exponentNum = 2.8;
M_nn30nd28.EarlyVis_spec.divNorm_exponentDen = 3.0;
M_nn30nd28.EarlyVis_spec.divNorm_rateScale_sps = 10;
M_nn30nd28 = DMPL_prepareSpecs(M_nn30nd28);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;
neuron_orient_idx = FindClosestIdx(M_nn30nd28.EarlyVis_spec.domain_orient_deg,neuron_orient_deg);
neuron_spFreq_idx = FindClosestIdx(M_nn30nd28.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);

maxSize_nd215 = max(abs(M_nn30nd28.stim_spec.gridX_deg(:)))*2;


%% Size tuning
specs_Size.contrast_pcent = 100;
specs_Size.neuron_orient_deg = neuron_orient_deg;
specs_Size.neuron_spFreq_cpd = neuron_spFreq_cpd;
data_Size = Explore_Size(M_nn30nd28,specs_Size);

% data_Size.plot('complex');
mCRF = data_Size.peak_x(complex);
fprintf('measured CRF = %f\n',mCRF);


%% Contrast sensitivity with Surround-suppression
specs_SurrCon.neuron_orient_deg = neuron_orient_deg;
specs_SurrCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
specs_SurrCon.contrast_pcent = 100;
specs_SurrCon.center_diameter_deg   = mCRF;
specs_SurrCon.surround_outer_diameter_out_deg = maxSize_nd215;
specs_SurrCon.surround_inner_diameter_out_deg = mCRF;
data_SurrCon = Explore_SurrContrast_on_Contrast(M_nn30nd28,specs_SurrCon);

data_SurrCon.plot('complex');

