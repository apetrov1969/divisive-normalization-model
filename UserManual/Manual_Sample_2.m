%% Sample script 2 for the Divisive Normalization model
%  Responses of the DNM to multiple visual stimuli from image files ('SampleImg_2a/b/c/d.bmp')

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt

%% Clean up for debugging
clear; close all; clc;


%% 1) Load image files
stimuli = zeros(64,64,4); % 64pixels x 64pixels == default image size of the DNM (see Table 2)

% A: Orientation = 120 deg, Frequency = 1.5 oct (2.83 cpd), Phase = 0
% B: Orientation = 135 deg, Frequency = 1.5 oct (2.83 cpd), Phase = 0
% C: Orientation = 135 deg, Frequency = 1.5 oct (2.83 cpd), Phase = 90
% C: Orientation = 135 deg, Frequency = 0.5 oct (1.41 cpd), Phase = 90

imageFileAddresses = {  'SampleImages/SampleImg_2a.bmp',...
                        'SampleImages/SampleImg_2b.bmp',...
                        'SampleImages/SampleImg_2c.bmp',...
                        'SampleImages/SampleImg_2d.bmp' };

for i=1:4
    loadedImg = double(imread(imageFileAddresses{i})); % 64x64
    if ndims(loadedImg)==3; loadedImg = mean(loadedImg,3); end; % from RGB to Gray-scale
    stimuli(:,:,i) = loadedImg;
end

stimuli = stimuli./255; % Normalize values between 0 and 1


%% 2) Set-up a "model" structure M with the standard parameter set (see Table 2)
M = DMPL_prepareSpecs(DMPL_defaultSpecs([])); % Set-up the DNM based on the standard parameter set (see Table 2)


%% 3) Apply the DNM to the stimuli
resp_DivNorm = DMPL_EarlyVis_FiringRate(M, stimuli); % Responses of 300 (12x5x1x5) different neurons to 2 stimuli:
  % [Orientation, Frequency, RF-location, Phase, Stimulus] 
  % 12 orientations:        M.EarlyVis_spec.domain_orient_deg   (0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165)
  %  5 spatial-frequencies: M.EarlyVis_spec.domain_spFreq_l2cpd (0, 0.5, 1, 1.5, 2)
  %  1 RF-locations:        M.EarlyVis_spec.domain_rfLoc_deg    ([0 0])
  %  5 phases: complex and Simple x 4 phases (0, 90, 180, 270)
  %  4 stimuli: 'SampleImg_2a/b/c/d.bmp'


%% 4) Response of the DNM to the stimuli
resultsOri = [  squeeze(resp_DivNorm(:,4,1,1,1))... [12 orientations]
                squeeze(resp_DivNorm(:,4,1,1,2))... [12 orientations]
                squeeze(resp_DivNorm(:,4,1,1,3))... [12 orientations]
                squeeze(resp_DivNorm(:,2,1,1,4))... [12 orientations]
             ]; % [12 orientations x 4 stimuli]

resultsSpf = [  squeeze(resp_DivNorm( 9,:,1,1,1))'... [5 spatial-frequencies]
                squeeze(resp_DivNorm(10,:,1,1,2))'... [5 spatial-frequencies]
                squeeze(resp_DivNorm(10,:,1,1,3))'... [5 spatial-frequencies]
                squeeze(resp_DivNorm(10,:,1,1,4))'... [5 spatial-frequencies]
             ]; % [5 spatial-frequencies x 4 stimuli]

resultsPha = [  squeeze(resp_DivNorm( 9,4,1,:,1))... [5 neuron-types]
                squeeze(resp_DivNorm(10,4,1,:,2))... [5 neuron-types]
                squeeze(resp_DivNorm(10,4,1,:,3))... [5 neuron-types]
                squeeze(resp_DivNorm(10,2,1,:,4))... [5 neuron-types]
             ]; % [5 neuron-types x 4 stimuli]

% Plot graphs
figure('Name','Sample 2');
set(gcf,'pos',[100 100 1200 900]);

for i=1:4
    subplot(4,4,4*i-4+1,'align');
    imshow(stimuli(:,:,i));

    subplot(4,4,4*i-4+2,'align');
    bar(M.EarlyVis_spec.domain_orient_deg, resultsOri(:,i), 'b');
	xlabel('Tuned orientation (deg)'); ylabel('Firing rate (spikes/sec)');
    axis([-15,180,0,50]); set(gca,'xtick',[0 45 90 135]);

    subplot(4,4,4*i-4+3,'align');
    bar(M.EarlyVis_spec.domain_spFreq_l2cpd, resultsSpf(:,i), 'r');
    xlabel('Tuned frequency (oct)'); ylabel('Firing rate (spikes/sec)');
    axis([-0.5,2.5,0,50]);

    subplot(4,4,4*i-4+4,'align');
    bar(resultsPha(:,i), 'black');
    xlabel('Tuned phase (deg)'); ylabel('Firing rate (spikes/sec)');
    ylim([0,50]); set(gca,'xticklabel',{'Complex' '0' '90' '180' '270'});
end

