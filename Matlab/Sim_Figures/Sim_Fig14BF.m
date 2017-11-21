%% Simulation Experiments using the standardized Divisive Normalization model
%  Effect of luminance contrast on bandwidths in the orientation/spatial-frequency domains
%  m_{\Theta_w} = 40, h_{F_w} = 1.0

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
M_BwCont = DMPL_defaultSpecs([],0,model_level);
M_BwCont.EarlyVis_spec.fwhh_pool_surr_spFreq_oct = 1.0; %2.0
M_BwCont.EarlyVis_spec.fwmh_pool_surr_orient_deg = 40;  %60
M_BwCont = DMPL_prepareSpecs(M_BwCont);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;
neuron_orient_idx = FindClosestIdx(M_BwCont.EarlyVis_spec.domain_orient_deg,neuron_orient_deg);
neuron_spFreq_idx = FindClosestIdx(M_BwCont.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);

maxSize_BwCont = max(abs(M_BwCont.stim_spec.gridX_deg(:)))*2;


%% Orientation/Spatial-frequency bandwidths by contrast

    %% Orientation tuning by contrast
    specs_OriByCon.neuron_orient_deg = neuron_orient_deg;
    specs_OriByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_OriByCon.diameter_deg = maxSize_BwCont;
    data_OriByCon = Explore_Orient_by_Contrast(M_BwCont,specs_OriByCon);

    data_OriByCon.plot('complex');

    %% Spatial-frequency tuning by contrast
    specs_SpfByCon.neuron_orient_deg = neuron_orient_deg;
    specs_SpfByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpfByCon.diameter_deg = maxSize_BwCont;
    data_SpfByCon = Explore_SpFreq_by_Contrast(M_BwCont,specs_SpfByCon);

    data_SpfByCon.plot('complex');
