%% Simulation Experiments using the standardized Divisive Normalization model
%  A modified parameter set for size tuning function

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
M_Hole = DMPL_defaultSpecs([],0,model_level);
M_Hole.stim_spec.imageSize_pix = [128,128];

M_Hole.EarlyVis_spec.divNorm_exponentDen = 2.5; % 2.0
M_Hole.EarlyVis_spec.divNorm_baselineConst = 0.0050; % 0.0200
M_Hole.EarlyVis_spec.divNorm_semisaturConst= 0.0400; % 0.1000
M_Hole.EarlyVis_spec.divNorm_rateScale_sps = 20; % 40

M_Hole = DMPL_prepareSpecs(M_Hole);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;

maxSize_Hole = max(abs(M_Hole.stim_spec.gridX_deg(:)))*2;


%% Subsection: Size tuning

    %% Size tuning
    specs_Size.contrast_pcent = 100;
    specs_Size.neuron_orient_deg = neuron_orient_deg;
    specs_Size.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_Size = Explore_Size(M_Hole,specs_Size);

    data_Size.plot('complex');
    mCRF = data_Size.peak_x(complex);
    fprintf('measured CRF = %f\n',mCRF);

    %% Hole-size tuning for Revision
    specs_HoleSizeByConRev.neuron_orient_deg = neuron_orient_deg;
    specs_HoleSizeByConRev.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_HoleSizeByConRev.contrast_pcent = 2;
    specs_HoleSizeByConRev.outer_diameter_deg = 3;%maxSize_Hole;
    data_HoleSizeByConRev = Explore_SizeHole(M_Hole,specs_HoleSizeByConRev);

    data_HoleSizeByConRev.plot('complex');

