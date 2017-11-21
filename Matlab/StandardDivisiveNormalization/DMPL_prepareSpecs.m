function  M = DMPL_prepareSpecs(M, levelIn, levelOut)

%DMPL_prepareSpecs - Fill in the derivable fields of Dimple model specs
%
% M = DMPL_prepareSpecs(M)
% M = DMPL_prepareSpecs(M, levelIn, levelOut)
%
% The input M is a specification structure (aka "model") generated by
% DMPL_defaultSpecs.m or some equivalent function.
% DMPL_prepareSpecs fills in certain additional fields that are needed for
% the operation of the model but are not free parameters because they are
% entirely determined by other parameters.
%
% DMPL_prepareSpecs is organized as a pipeline whose levels follow the
% levels of the DMPL_defaultSpecs pipeline. 
% Type 'helpwin DMPL_defaultSpecs' for details and examples.
% Each stage of the DMPL_prepareSpecs pipeline is implemented as a separate
% function as follows:
%  level 0: Only a descriptor string with a timestamp and version info
%  level 1: See DMPL_stim_prepareSpecs
%  level 2: See DMPL_EarlyVis_prepareSpecs
%  level 3: Not currently implemented. Will be added later...
%
% Input arguments:  ----------------------
% M -- Spec structure generated by DMPL_defaultSpecs or equivalent
% levelIn -- Optional entry point on the pipeline. Default = 0
% levelOut -- Optional exit point. Default = M0.model_level
%
% Return value:  -------------------------
% M -- An updated structure with fields determied by
%   DMPL_stim_prepareSpecs.m, DMPL_EarlyVis_prepareSpecs.m, etc.
%
% Documentation is available in .../Dimple0/docs/Dimple_frontend_eqs.tex
%
% Example 1: Generate everything from scratch:
%  M = DMPL_defaultSpecs; M = DMPL_prepareSpecs(M)
%    M =        descriptor: 'DMPL_defaultSpecs ver. 0sa, 22-Jun-2013 17:46:43'
%              model_level: 2
%                stim_spec: [1x1 struct]
%            EarlyVis_spec: [1x1 struct]
%                     ....: ....   % more default specs added by higher levels
%                 preparer: 'DMPL_prepareSpecs, 22-Jun-2013 17:46:43'
%            prepare_level: 2
%         EarlyVis_filters: [1x1 struct]
%     EarlyVis_poolKernels: [1x1 struct]
%         EarlyVis_divNorm: [1x1 struct]
%     EarlyVis_calibration: [1x1 struct]
%                     ....: ....   % more prepared specs added by higher levels
% 
% Example 2: Change some parameter, then redo the *relevant* preparation
%  % For concreteness, suppose we want to sharpen the orientation bandwidth:
%  M.EarlyVis_spec.fwhh_bandwidth_orient_deg = 25;
%  % We must redo all preparation from level 2 upward:
%  M = DMPL_prepareSpecs(M,2)
%  % On the surface M looks identical as before, but important changes
%  % have been made inside M.EarlyVis_filters and M.EarlyVis_calibration.
%
% See also DMPL_defaultSpecs, DMPL_stim_prepareSpecs

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.1 2015-12-17 TS: Modified
% 1.0.0 2013-06-16 AAP: Wrote it from scratch


%% Supply default input args if needed
if (nargin<2); levelIn  = 0  ; end     % default
if (nargin<3); levelOut = M.model_level; end

%- Check level args for consistency
assert(levelIn <= levelOut)


%% Level 0:  Base class   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.preparer = sprintf('%s, %s', mfilename(), datestr(now));
if (levelIn <= 0)
    M.prepare_level = 0;   % incremented by each pipeline stage below
%else % do nothing here and continue with higher levels
end
if (levelOut <= 0); return; end   % skip the remaining levels


%% Level 1:  Stimulus specs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (levelIn <= 1)
    M = DMPL_stim_prepareSpecs(M);
%else % do nothing here and continue with higher levels
end
if (levelOut <= 1); return; end   % skip the remaining levels


%% Level 2:  EarlyVis specs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (levelIn <= 2)
    M = DMPL_EarlyVis_prepareSpecs(M);
%else % do nothing here and continue with higher levels
end
if (levelOut <= 2); return; end   % skip the remaining levels


%% Level 3:  To be added later ...


%% Return M
end  %%% of file
