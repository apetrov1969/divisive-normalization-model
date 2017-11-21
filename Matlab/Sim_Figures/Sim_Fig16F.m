%% Simulation Experiments using the standardized Divisive Normalization model
%  A modified parameter set for cross-orientation suppression

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
M_Strong = DMPL_defaultSpecs([],0,model_level);

M_Strong.EarlyVis_spec.fwmh_pool_surr_orient_deg = 90; %60
M_Strong.EarlyVis_spec.fwhh_pool_surr_spFreq_oct = 1.0; %3.0

M_Strong.EarlyVis_spec.divNorm_rateScale_sps = 50;

M_Strong = DMPL_prepareSpecs(M_Strong);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;

maxSize_Small = max(abs(M_Strong.stim_spec.gridX_deg(:)))*2;


%% Subsection: Spatial-frequency tuning

    %% Spatial-frequency tuning
%     specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
%     specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
%     specs_SpFreq.contrast_pcent = 50;
%     specs_SpFreq.diameter_deg = maxSize_Small;
%     data_SpFreq = temp_Explore_SpFreq(M_Small,specs_SpFreq);
% 
%     % data_SpFreq.plot('complex');
%     
%     idx_neuron_spFreq = FindClosestIdx(M_Small.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);
%     peak_spf_oct = data_SpFreq.marking_x(complex,1,idx_neuron_spFreq);
%     peak_spf_cpd = 2^peak_spf_oct %#ok
    
    peak_spf_cpd = 2;
    
    
%% Subsection: X-orientation-suppression
    

    %% Orientation tuning of X-suppression
%     specs_XOri.neuron_orient_deg = neuron_orient_deg;
%     specs_XOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
%     specs_XOri.base_spFreq_cpd = peak_spf_cpd;
%     specs_XOri.contrast1_pcent = 15;
%     specs_XOri.contrast2_pcent = 25;
%     specs_XOri.diameter_deg = maxSize_Small;
%     data_XOri = Explore_XOrient(M_Small,specs_XOri);
% 
%     data_XOri.plot('complex');

    %% Spatial-frequency of X-suppression
    specs_XSpf.neuron_orient_deg = neuron_orient_deg;
    specs_XSpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_XSpf.base_spFreq_cpd = peak_spf_cpd;
    specs_XSpf.contrast1_pcent = 10;
    specs_XSpf.contrast2_pcent = 25;
    specs_XSpf.diameter_deg = maxSize_Small;
    data_XSpf = Explore_XSpFreq(M_Strong,specs_XSpf);

    data_XSpf.plot('complex');
    

    
    