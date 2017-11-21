function  M = DMPL_stim_defaultSpecs(M)
%DMPL_stim_defaultSpecs - Default stimulus specs of the Dimple model
%
% M1 = DMPL_stim_defaultSpecs(M0)
%
% DMPL_stim_defaultSpecs generates default specifications for the images
% passed as inputs to the DMPL model.
%
% DMPL_stim_defaultSpecs implements level 1 of the DMPL_defaultSpecs pipeline.
%
% Input arguments:  ----------------------
% M0 -- Seed struct produced by (or compatible with) DMPL_defaultSpecs level 0.
%
% Return value:  -------------------------
% M1 -- A copy of M0 plus one additional (or modified) field:
%  M.stim_spec -- Itself a struct with fields shown in the example below.
%
% Example:
% M0 = DMPL_defaultSpecs([],0,0);   % make a level-0 seed structure
% M1 = DMPL_stim_defaultSpecs(M0) , stim_spec = M1.stim_spec
%   M1 = 
%     descriptor: 'DMPL_defaultSpecs ver. 0sa, 15-Jun-2013 17:58:20'
%    model_level: 1            % incremented from M0.model_level
%      stim_spec: [1x1 struct]
%
%   stim_spec = 
%           descriptor: 'DMPL_stim_defaultSpecs, 15-Jun-2013 18:12:03'
%     coordinatesStyle: 'xy'
%        imageSize_pix: [64 64]
%          degPerPixel: 0.0450
%               minLum: 0
%                bgLum: 0.5000
%               maxLum: 1
%
% See also DMPL_defaultSpecs, DMPL_stim_prepareSpecs, DMPL_EarlyVis_defaultSpecs

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.2 2015-12-17 TS: Modified
% 1.0.1 2013-06-17 AAP: Envorce that (bgLum==(maxLum-minLum)/2) is always true.
% 1.0.0 2013-06-16 AAP: Wrote it based on Dimple0k/DMPL_EarlyVis_Specs.m


%% Check seed struct for consistency
%assert(isstruct(M))
assert(M.model_level >= 0)


%% Default stimulus specs
stim_spec.descriptor = sprintf('%s, %s', mfilename(), datestr(now));

stim_spec.coordinatesStyle = 'ij';
%'ij'(i:top to bottom, j:left to right), 0deg='|', 45deg='/'
%'xy'(y:bottom to top, x:left to right), 0deg='|', 45deg='\'
% The i-axis is corresponded with the y-axis with reflection
% The j-axis is corresponded with the x-axis
% e.g. imagesc(I); % Show the image I.
%      size(I);    % Show ranges of the i- and j-axes
%      axis('xy'); % This command reflect the image upside-down.

stim_spec.imageSize_pix = [64,64];

stim_spec.degPerPixel = 0.0450;  % proxy for viewing distance, etc.


%% Luminance-related specs
% The background luminance bgLum MUST be halfway between minLum and maxLum.
% This is a theoretical commitment of the model.  The image isn't really
% defined in terms of luminance but in terms of Michelson contrast.
% The model assumes that light adaptation, pupil dilation/contraction, and
% other gain-control mechanisms have taken place, such that bgLum reflects
% the prevailing average illumination level in the environment.
%
% If you want to model an experiment that violates the equality:
%   bgLum - minLum == maxLum - bgLum,
% then you must adjust minLum and/or maxLum so that your halfway point
% coincides with the background luminance of your stimuli.

stim_spec.minLum = 0;   % minimum luminance for a *set* of stimuli
stim_spec.maxLum = 1;   % maximum luminance for a *set* of stimuli
stim_spec.bgLum  = (stim_spec.maxLum - stim_spec.minLum)/2;


%% Add to model specs and update model level
M.stim_spec = stim_spec;
M.model_level = max(M.model_level,1);


%%% Return M
end  %%% of file
