function  M = DMPL_EarlyVis_defaultSpecs(M)
%DMPL_EarlyVis_defaultSpecs - Default stimulus specs of the Dimple model
%
% M1 = DMPL_EarlyVis_defaultSpecs(M0)
%
% DMPL_EarlyVis_defaultSpecs generates default specifications for the
% virtual neurons.
%
% See also DMPL_defaultSpecs, DMPL_stim_prepareSpecs, DMPL_EarlyVis_defaultSpecs

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.1 2015-12-17 TS: Modified


%% Check seed struct for consistency
assert(isstruct(M))
assert(isfield(M,'stim_spec'))  % produced by DMPL_stim_defaultSpecs at level 1


%% Default EarlyVis specs
EarlyVis_spec.descriptor = sprintf('%s, %s', mfilename(), datestr(now));


%% The types of the neurons
EarlyVis_spec.num_cellType = 5; % 1 Complex cell + 4 Simple cells


%% The orientation preferences of the stimulus-drives
num_orient = 12;  % number of orientation channels
EarlyVis_spec.num_orient = num_orient;
EarlyVis_spec.step_orient_deg = 180 / num_orient;

domain_orient_deg = 0:EarlyVis_spec.step_orient_deg:179.999;
assert(length(domain_orient_deg)==num_orient);
EarlyVis_spec.domain_orient_deg = domain_orient_deg;

%M.stim_spec.coordinatesStyle = 'ij':
%    i-axis: top to bottom
%    j-axis: left to right
%    Rotation: Clockwise
%    0deg='||', 45deg='//', 90='=', 135='\\'
%M.stim_spec.coordinatesStyle = 'xy':
%    y-axis:bottom to top
%    x-axis:left to right
%    Rotation: Counter-clockwise
%    0deg='||', 45deg='\\', 90='=', 135='//'


%% The spatial-frequency ("spFreq") preferences of the stimulus-drives
% f cycles/deg = log2(f) octaves
% o octaves = 2^o cycles/deg
num_spFreq = 5;
num_channel_spFreq = num_spFreq+2; % additional very low/high frequency channels
EarlyVis_spec.num_spFreq = num_spFreq;
EarlyVis_spec.num_channel_spFreq = num_channel_spFreq;
EarlyVis_spec.step_spFreq_oct = 0.5;

domain_channel_spFreq_l2cpd = (-1:num_channel_spFreq-2) .* EarlyVis_spec.step_spFreq_oct;
assert(length(domain_channel_spFreq_l2cpd)==num_channel_spFreq);

EarlyVis_spec.domain_channel_spFreq_l2cpd =    domain_channel_spFreq_l2cpd; % in log2 cycles per degree
EarlyVis_spec.domain_channel_spFreq_cpd   = 2.^domain_channel_spFreq_l2cpd; % in cycles per degree

EarlyVis_spec.domain_spFreq_l2cpd =    domain_channel_spFreq_l2cpd(2:num_channel_spFreq-1); % in log2 cycles per degree
EarlyVis_spec.domain_spFreq_cpd   = 2.^domain_channel_spFreq_l2cpd(2:num_channel_spFreq-1); % in cycles per degree


%% The orientation bandwidths (full width at half height) of the stimulus drives in degree
EarlyVis_spec.fwhh_bandwidth_orient_deg = 40.0;


%% The spatial-frequency bandwidths (full width at half height) of the stimulus-drives in octave
EarlyVis_spec.fwhh_bandwidth_spFreq_oct = 1.5;


%% The type of the image filters applied to the input images for channels
% 'gabor': Gabor filter
EarlyVis_spec.filters_type = 'gabor';


%% Phase invariance method for the suppressive-drive
% 'quadrature': Both cos and sin at individual pixels (Mathematically, better)
% 'cos-only': Cos only (Mathematically, approximation. Computationally, faster)
% 'sin-only': Sin only (Mathematically, approximation. Computationally, faster)
EarlyVis_spec.phaseInvariance_method = 'quadrature';


%% The xy-spatial positions of the neurons
% EarlyVis_spec. domain_rfLoc_deg(1) must be at the center [0 0]
EarlyVis_spec.num_rfLoc = 1;
EarlyVis_spec. domain_rfLoc_deg = [0 0]; % [y,x], locations of cells' receptive fields
assert(EarlyVis_spec.num_rfLoc == size(EarlyVis_spec.domain_rfLoc_deg,1));


%% The kernels of the xy-spatial pooling for the suppressive-drives
%- The sizes of the xy-spatial pooling for the suppressive-drives
%   The surround xy-spatial pooling for the suppressive-drives
%   at the preferred spatial frequency of M.EarlyVis_spec.calibr_spFreq_idx
EarlyVis_spec.fwhh_pool_surr_xy_cyc = 2.0; 


%- The sizes of the xy-spatial pooling <-> The prefeered spatial-frequencies
%  'constant': Constant pooling sizes for all spatial-frequencies
%  'linear': Linearly-correlated pooling sizes with spf
EarlyVis_spec.pool_xy_method = 'linear'; 


%% The kernels for integrating the channels across orientations for the suppressive-drives
% Full width at middle height of von Mises distribution in degree (0~90)
% See VonMisesDistrib_FWMH_to_Kappa.m
EarlyVis_spec.fwmh_pool_surr_orient_deg = 60;


%% The kernels for integrating the channels across spatial-frequencies for the suppressive-drives
% Full width at half height in octave
EarlyVis_spec.fwhh_pool_surr_spFreq_oct = 2.0;
EarlyVis_spec.shift_pool_surr_spFreq_oct = 0.0;


%% Parameters of the divisive normalization equation
EarlyVis_spec.divNorm_exponentNum = 2.0; % for the stimulus-drives in the numerator
EarlyVis_spec.divNorm_exponentDen = 2.0; % for the suppressive-drives in the denominator

EarlyVis_spec.divNorm_baselineConst  = 0.02; % constant in the numerator
EarlyVis_spec.divNorm_semisaturConst = 0.1;  % constant in the denominator

EarlyVis_spec.divNorm_rateScale_sps = 40;  % for scaling the output of the model (firing-rate)


%% Calibration related specs
%- Orientation for the calibration: vertical ('||', 0deg)
cal_ori_deg = 0; % deg
cal_idx_ori = find(EarlyVis_spec.domain_orient_deg==cal_ori_deg);
assert(~isempty(cal_idx_ori));
EarlyVis_calibration.calibr_orient_deg = cal_ori_deg;
EarlyVis_calibration.calibr_orient_chan_idx = cal_idx_ori;
EarlyVis_calibration.calibr_orient_cell_idx = cal_idx_ori;

%- Frequency for the calibration: 2 cycles/deg (= 1 octaves)
cal_spf_cpd = 2.0;
cal_chan_idx_spf = find(EarlyVis_spec.domain_channel_spFreq_cpd==cal_spf_cpd);
assert(~isempty(cal_chan_idx_spf));
cal_cell_idx_spf = find(EarlyVis_spec.domain_spFreq_cpd==cal_spf_cpd);
assert(~isempty(cal_cell_idx_spf));
EarlyVis_calibration.calibr_spFreq_cpd  = cal_spf_cpd;
EarlyVis_calibration.calibr_spFreq_chan_idx = cal_chan_idx_spf;
EarlyVis_calibration.calibr_spFreq_cell_idx = cal_cell_idx_spf;

%- Position for the calibration: center (assumed to have index 1)
EarlyVis_calibration.calibr_rfLoc_cell_idx  = 1;

%- Set temporal calibration factors
EarlyVis_calibration.calibrFactor_numerator_simp = 1;
EarlyVis_calibration.calibrFactor_numerator_comp = 1;
EarlyVis_calibration.calibrFactor_denominator = 1;


%% Noise related specs
% To be added later


%% Add to model specs and update model level
M.EarlyVis_spec = EarlyVis_spec;
M.EarlyVis_calibration = EarlyVis_calibration;
M.model_level = max(M.model_level,2);


%%% Return M
end  %%% of file



%% Removed
%- The sizes of the xy-spatial pooling for the stimulus-drives
%   The center xy-spatial pooling for the stimulus-drives
%   at the preferred spatial frequency of M.EarlyVis_spec.calibr_spFreq_idx
%EarlyVis_spec.fwhh_pool_cent_xy_deg = 2.0;

