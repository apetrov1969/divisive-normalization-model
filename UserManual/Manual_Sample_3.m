%% Sample script 3 for the Divisive Normalization model
%  Responses of the DNM with modified parameters to a visual stimulus

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt

%% Clean up for debugging
clear; close all; clc;


%% 1) Load an image file
loadedImg = double(imread('SampleImages/SampleImg_3.bmp')); % 256x256
if ndims(loadedImg)==3; loadedImg = mean(loadedImg,3); end; % from RGB to Gray-scale
stimulus = loadedImg./255; % Normalize values between 0 and 1


%% 2) Set-up a "model" structure M with modified parameters
M = DMPL_defaultSpecs([]); % The standard parameter set of the DNM (see Table 2)

% Modify stimulus parameters
M.stim_spec.imageSize_pix = [size(stimulus,1),size(stimulus,2)]; % Image size
M.stim_spec.degPerPixel = 0.02; % Image resolution

% Modify the DNM parameters
M.EarlyVis_spec.divNorm_exponentNum = 2.5; % $n_n$
M.EarlyVis_spec.divNorm_exponentDen = 2.5; % $n_d$

M.EarlyVis_spec.fwhh_bandwidth_orient_deg = 50; % $h_{\theta}$
M.EarlyVis_spec.fwhh_bandwidth_spFreq_oct = 2; % $h_f$
  
M.EarlyVis_spec.fwmh_pool_surr_orient_deg = 90; % $h_{\Theta}$
M.EarlyVis_spec.fwhh_pool_surr_spFreq_oct = 1; % $h_F$
M.EarlyVis_spec.shift_pool_surr_spFreq_oct = -1; % $\myu_F$

M.EarlyVis_spec.divNorm_baselineConst = 0.0; % $\beta$
M.EarlyVis_spec.divNorm_semisaturConst= 0.05; % $\alpha$

M.EarlyVis_spec.divNorm_rateScale_sps = 4; % $R$

M.EarlyVis_spec.domain_rfLoc_deg = [... 
                                    1.2*sin(1*pi/6) 1.2*cos(1*pi/6);... [y x]
                                    1.2*sin(2*pi/6) 1.2*cos(2*pi/6);... [y x]
                                    1.2*sin(3*pi/6) 1.2*cos(3*pi/6);... [y x]
                                    1.2*sin(4*pi/6) 1.2*cos(4*pi/6);... [y x]
                                    1.2*sin(5*pi/6) 1.2*cos(5*pi/6);... [y x]
                                    ...
                                    0 0.6;... [y x]
                                    0 1.0;... [y x]
                                    0 1.4;... [y x]
                                    0 1.8;... [y x]
                                    0 2.2;... [y x]
                                    ...
                                    -1.07 0;... [y x]
                                    -1.22 0;... [y x]
                                    -1.33 0;... [y x]
                                    -1.43 0;... [y x]
                                   ]; % Receptive field centers

% Set-up the model structure M based on the modified parameters
M = DMPL_prepareSpecs(M);


%% 3) Apply the DNM to the stimulus
resp_DivNorm = DMPL_EarlyVis_FiringRate(M, stimulus); % Responses of 300 (12x5x1x5) different neurons to 2 stimuli:
  % [Orientation, Frequency, RF-location, Phase, Stimulus] 
  % 12 orientations:        M.EarlyVis_spec.domain_orient_deg   (0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165)
  %  5 spatial-frequencies: M.EarlyVis_spec.domain_spFreq_l2cpd (0, 0.5, 1, 1.5, 2)
  % 14 RF-locations:        M.EarlyVis_spec.domain_rfLoc_deg
  %  5 phases: complex and Simple x 4 phases (0, 90, 180, 270)
  %  1 stimulus (This dimension is omitted)


%% 4) Response of the DNM to the stimulus
resultsOri = squeeze(resp_DivNorm(:,3, 1: 5,1)); % [12 orientations, 1:5 locations]
resultsSpf = squeeze(resp_DivNorm(1,:, 6:10,1)); % [ 5 spatial-frequencies, 6:10 locations]
resultsPha = squeeze(resp_DivNorm(7,3,11:14,:)); % [ 5 neuron-types, 11:15 locations]
    
% Plot graphs
figure('Name','Sample 3');
set(gcf,'pos',[100 100 1200 900]);

for i=1:5
    subplot(3,5,i,'align');
    bar(M.EarlyVis_spec.domain_orient_deg, resultsOri(:,i), 'b');
 	xlabel('Tuned orientation (deg)'); ylabel('Firing rate (spikes/sec)');
    axis([-15,180,0,30]); set(gca,'xtick',[0 45 90 135]);
    
    subplot(3,5,i+5,'align');
    bar(M.EarlyVis_spec.domain_spFreq_l2cpd, resultsSpf(:,i), 'r');
    xlabel('Tuned frequency (oct)'); ylabel('Firing rate (spikes/sec)');
    axis([-0.5,2.5,0,30]);
    
end

for i=1:4
    subplot(3,5,i+10,'align');
    bar(resultsPha(i,:), 'black');
    xlabel('Tuned phase (deg)'); ylabel('Firing rate (spikes/sec)');
    ylim([0,50]); set(gca,'xticklabel',{'Complex' '0' '90' '180' '270'});
end

% Show receptive fieled centers on the stimulus
DMPL_Supplement_ShowLocationsInImage(M, stimulus);