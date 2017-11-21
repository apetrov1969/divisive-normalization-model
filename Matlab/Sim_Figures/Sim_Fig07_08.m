%% Simulation Experiments using the standardized Divisive Normalization model
%  Stimuli with higher resolution (0.02 deg/pixel) for Size tuning

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
% The regular parameter set (for finding sub-optimal orientations/spatial-frequencies)
model_level = 2;
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small = DMPL_prepareSpecs(M_Small);

% Higher resolution
model_level = 2;
M_HighReso = M_Small;
M_HighReso.stim_spec.imageSize_pix = [256,256];
M_HighReso.stim_spec.degPerPixel = 0.02;
M_HighReso = DMPL_prepareSpecs(M_HighReso);


%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;
neuron_orient_idx = FindClosestIdx(M_HighReso.EarlyVis_spec.domain_orient_deg,neuron_orient_deg);
neuron_spFreq_idx = FindClosestIdx(M_HighReso.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);

maxSize_Small = max(abs(M_Small.stim_spec.gridX_deg(:)))*2;
maxSize_HighReso = max(abs(M_HighReso.stim_spec.gridX_deg(:)))*2;


%% Subsection: Size tuning

    %% Size tuning
    specs_Size.contrast_pcent = 100;
    specs_Size.neuron_orient_deg = neuron_orient_deg;
    specs_Size.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_Size = Explore_Size(M_Small,specs_Size);

    data_Size.plot('complex');
    mCRF = data_Size.peak_x(complex);
    fprintf('measured CRF = %f\n',mCRF);

    %% Size tuning (orientation)
        % Finding sub-optimal orientation
        subopt_specs_Orient.contrast_pcent = 100;
        subopt_specs_Orient.diameter_deg = mCRF;
        subopt_specs_Orient.neuron_orient_deg = neuron_orient_deg;
        subopt_specs_Orient.neuron_spFreq_cpd = neuron_spFreq_cpd;
        subopt_data_Orient = Explore_Orient(M_Small,subopt_specs_Orient);

    specs_SizeByOri.neuron_orient_deg = neuron_orient_deg;
    specs_SizeByOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SizeByOri.fwhh_orient_deg = subopt_data_Orient.bandwidth_deg(complex);
    data_SizByOri = Explore_Size_by_Orient(M_HighReso,specs_SizeByOri);

    data_SizByOri.plot('complex');

    %% Size tuning (spatial-frequency)
        % Finding sub-optimal spatial-frequency
        subopt_specs_SpFreq.contrast_pcent = 100;
        subopt_specs_SpFreq.diameter_deg = mCRF;
        subopt_specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
        subopt_specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
        subopt_data_SpFreq = Explore_SpFreq(M_Small,subopt_specs_SpFreq);

    specs_SizeBySpf.neuron_orient_deg = neuron_orient_deg;
    specs_SizeBySpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SizeBySpf.neuron_spFreq_idx = neuron_spFreq_idx;
    specs_SizeBySpf.spFreq_low_oct  = subopt_data_SpFreq.marking_x(complex,2,specs_SizeBySpf.neuron_spFreq_idx);
    specs_SizeBySpf.spFreq_high_oct = subopt_data_SpFreq.marking_x(complex,3,specs_SizeBySpf.neuron_spFreq_idx);
    data_SizBySpf = Explore_Size_by_SpFreq(M_HighReso,specs_SizeBySpf);

    data_SizBySpf.plot('complex');