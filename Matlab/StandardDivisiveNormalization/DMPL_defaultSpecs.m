function  M = DMPL_defaultSpecs(M, levelIn, levelOut)

%DMPL_defaultSpecs - Default specifications of the Dimple model
%
% M = DMPL_defaultSpecs()
% M = DMPL_defaultSpecs(M0, levelIn, levelOut)
%
% DMPL_defaultSpecs generates default specifications for the Dual-Pathway
% Model of Perceptual Learning (Dimple).  The spec structure M is
% constructed incrementally, in a pipeline divided into several 'levels'
% as follows:
%  level 0: Only a descriptor string with a timestamp and version info
%  level 1: See DMPL_stim_defaultSpecs
%  level 2: See DMPL_EarlyVis_defaultSpecs
%  level 3: Not currently implemented. Will be added later...
%
% Input arguments:  ----------------------
% M0 -- Optional seed structure.  When not supplied (or []), default
%   specifications are generated from scratch.
%
% levelIn -- Optional entry point on the pipeline. Default = 0
% levelOut -- Optional exit point. Default = Inf
%
% Return value:  -------------------------
% M -- A specification (or "spec") structure with fields determied by
%   DMPL_stim_defaultSpecs.m, DMPL_EarlyVis_defaultSpecs.m, etc.
%
% Documentation is available in .../Dimple0/docs/Dimple_frontend_eqs.tex
%
% Example 1: Typical usage -- generate everything from scratch:
% M = DMPL_defaultSpecs
%   M = 
%        descriptor: 'DMPL_defaultSpecs ver. 0sa, 15-Jun-2013 18:13:32'
%       model_level: 2
%         stim_spec: [1x1 struct]     % level 1
%     EarlyVis_spec: [1x1 struct]     % level 2
%              ....: ....             % more specs added by higher levels
% 
% Example 2: Generate stimulus specs only and exit at level 1
% M1 = DMPL_defaultSpecs([],0,1)
%   M1 = 
%     descriptor: 'DMPL_defaultSpecs ver. 0sa, 15-Jun-2013 18:13:32'
%    model_level: 1
%      stim_spec: [1x1 struct]
%
% Example 3: Build on the existing level-1 specs in M1
% M = DMPL_defaultSpecs(M1,2)
%   % Produces the same output as Example 1 above
%
% See also DMPL_EarlyVis_defaultSpecs, DMPL_prepareSpecs

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2015-12-17 TS: Modified
% 1.0.0 2013-06-15 AAP: Wrote it from scratch


version = '0sa';

if (nargin<2)
    levelIn = 0;
    levelOut = 2;
elseif (nargin<3)
    levelOut = 2;
end


%% Supply default input args if needed
if (nargin<1 || isempty(M))    % no seed stucture is supplied?
    M = struct('descriptor',sprintf('%s ver. %s, %s', mfilename(), version, datestr(now)));
end

if (nargin<2); levelIn  = 0  ; end     % default
if (nargin<3); levelOut = Inf; end

%- Check level args for consistency
assert(levelIn <= levelOut)


%% Level 0:  Base class   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.model_level = 0;   % incremented by each pipeline stage below
if (levelOut <= 0); return; end   % skip the remaining levels


%% Level 1:  Stimulus specs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (levelIn <= 1)
    M = DMPL_stim_defaultSpecs(M);
%else % do nothing here and continue with higher levels
end
if (levelOut <= 1); return; end   % skip the remaining levels


%% Level 2:  EarlyVis specs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (levelIn <= 2)
    M = DMPL_EarlyVis_defaultSpecs(M);
%else % do nothing here and continue with higher levels
end
if (levelOut <= 2); return; end   % skip the remaining levels


%% Level 3:  TBC


%% Return M
end  %%% of file
