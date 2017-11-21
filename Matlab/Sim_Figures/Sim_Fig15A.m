%% Simulation Experiments using the standardized Divisive Normalization model
%  Comparison between sinusoidal/square gratings
%  h_f = 0.8, h_{F_w} = 0.4

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
M_SinSquare1 = DMPL_defaultSpecs([],0,model_level);
M_SinSquare1.stim_spec.imageSize_pix = [128,128]; %[64,64]

% Resolution of visual stimuli was increased for this simulation
% for clear distinction between sinusoidal/square gratings.
M_SinSquare1.stim_spec.degPerPixel = 0.02; %0.0450
M_SinSquare1 = DMPL_prepareSpecs(M_SinSquare1);

M_SinSquare2 = M_SinSquare1;
M_SinSquare2.EarlyVis_spec.fwhh_bandwidth_spFreq_oct = 0.8; %1.5;
M_SinSquare2.EarlyVis_spec.fwhh_pool_surr_spFreq_oct = 0.4; %2
M_SinSquare2 = DMPL_prepareSpecs(M_SinSquare2);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;
neuron_orient_idx = FindClosestIdx(M_SinSquare1.EarlyVis_spec.domain_orient_deg,neuron_orient_deg);
neuron_spFreq_idx = FindClosestIdx(M_SinSquare1.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);

maxSize_SinSquare = max(abs(M_SinSquare1.stim_spec.gridX_deg(:)))*2;


%% Subsection: M_SinSquare1

    %% Contrast sensitivity with sinusoidal/square waves
    specs_Contrast.neuron_orient_deg = neuron_orient_deg;
    specs_Contrast.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Contrast.diameter_deg = maxSize_SinSquare;
    data_Contrast1 = Explore_Contrast_by_SinSquare(M_SinSquare1,specs_Contrast);

    data_Contrast1.plot('complex');

    %% Spatial-frequency tuning with sinusoidal/square waves
    specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
    specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpFreq.contrast_pcent = 100;
    specs_SpFreq.diameter_deg = maxSize_SinSquare;
    data_SpFreq1 = Explore_SpFreq_by_SinSquare(M_SinSquare1,specs_SpFreq);

    data_SpFreq1.plot('complex');
    

%% Subsection: M_SinSquare2

    %% Contrast sensitivity with sinusoidal/square waves
    specs_Contrast.neuron_orient_deg = neuron_orient_deg;
    specs_Contrast.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Contrast.diameter_deg = maxSize_SinSquare;
    data_Contrast2 = Explore_Contrast_by_SinSquare(M_SinSquare2,specs_Contrast);

    data_Contrast2.plot('complex');

    %% Spatial-frequency tuning with sinusoidal/square waves
    specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
    specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpFreq.contrast_pcent = 100;
    specs_SpFreq.diameter_deg = maxSize_SinSquare;
    data_SpFreq2 = Explore_SpFreq_by_SinSquare(M_SinSquare2,specs_SpFreq);

    data_SpFreq2.plot('complex');