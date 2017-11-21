%% Sample script 1 for the Divisive Normalization model
%  Responses of the DNM to a visual stimulus from image fils ('SampleImg_1.bmp')

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt

%% Clean up for debugging
clear; close all; clc;


%% 1) Load an image file
stimulus = zeros(64,64); % 64pixels x 64pixels == default image size of the DNM (see Table 2)

%  Orientation = 15 deg, Frequency = 1 oct (2 cpd), Phase = 0
loadedImg = double(imread('SampleImages/SampleImg_1.bmp')); % 64x64
if ndims(loadedImg)==3; loadedImg = mean(loadedImg,3); end; % from RGB to Gray-scale
stimulus(:,:) = loadedImg./255; % Normalize values between 0 and 1


%% 2) Set-up a "model" structure M with the standard parameter set
M = DMPL_prepareSpecs(DMPL_defaultSpecs([])); % Set-up the DNM based on the standard parameter set (see Table 2)


%% 3) Apply the DNM to the stimulus
resp_DivNorm = DMPL_EarlyVis_FiringRate(M, stimulus); % Responses of 300 (12x5x1x5) different neurons to a stimulus:
  % [Orientation, Frequency, RF-location, Phase, Stimulus] 
  % 12 orientations:        M.EarlyVis_spec.domain_orient_deg   (0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165)
  %  5 spatial-frequencies: M.EarlyVis_spec.domain_spFreq_l2cpd (0, 0.5, 1, 1.5, 2)
  %  1 RF-locations:        M.EarlyVis_spec.domain_rfLoc_deg    ([0 0])
  %  5 phases: complex and Simple x 4 phases (0, 90, 180, 270)
  %  1 stimulus (This dimension is omitted)


%% 4) Response of the DNM to the stimulus
iOri = 2; % M.EarlyVis_spec.domain_orient_deg(2) = 15 deg
iSpf = 3; % M.EarlyVis_spec.domain_spFreq_l2cpd(iSpf) = 1 oct = 2 cpd
iRfl = 1; % M.EarlyVis_spec.domain_rfLoc_deg(iRfl) = [0 0]
iPhs = 1; % Complex cell
DMPL_Supplement_ShowFilter(M,iOri,iSpf);

resultsOri = squeeze(resp_DivNorm(:,iSpf,iRfl,iPhs)); % [12 orientations]
resultsSpf = squeeze(resp_DivNorm(iOri,:,iRfl,iPhs)); % [5 spatial-frequencies] 
resultsPha = squeeze(resp_DivNorm(iOri,iSpf,iRfl,:)); % [5 neuron-types, 2 stimulus]

% Plot graphs
figure('Name','Sample 1');
set(gcf,'pos',[100 100 700 700]);

subplot(2,2,1,'align');
imshow(stimulus);

subplot(2,2,2,'align');
bar(M.EarlyVis_spec.domain_orient_deg, resultsOri, 'b');
xlabel('Tuned orientation (deg)'); ylabel('Firing rate (spikes/sec)');
axis([-15,180,0,50]); set(gca,'xtick',[0 45 90 135]);

subplot(2,2,3,'align');
bar(M.EarlyVis_spec.domain_spFreq_l2cpd, resultsSpf, 'r');
xlabel('Tuned frequency (oct)'); ylabel('Firing rate (spikes/sec)');
axis([-0.5,2.5,0,50]);

subplot(2,2,4,'align');
bar(resultsPha, 'black');
xlabel('Tuned phase (deg)'); ylabel('Firing rate (spikes/sec)');
ylim([0,50]); set(gca,'xticklabel',{'Complex' '0' '90' '180' '270'});


