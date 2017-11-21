function linFWHH = Log_to_Linear_FWHH(logFWHH,center,base)

% Log_to_Linear_FWHH -- Convert log FWHH to linear FWHH bandwidth parameter
%
% linFWHH = Log_to_Linear_FWHH(logFWHH,center,base)
%
% Often full-width at half-height (FWHH) bandwidths of spatial-frequency filters
% are specified in octaves, which are on a logarithmic (base-2) scale.
% Log_to_Linear_FWHH converts to an equivalent FWHH on a linear scale.
% The linear bandwidth also depends on the center frequency as described below.
%
% Input arguments:  ----------------------
%  logFWHH -- Full-width at half-height in log units (typically octaves, i.e. base=2)
%  center  -- (A matrix of) the center(s) of the tuning curve(s).
%  base    -- Optional base of the logarithm.  Default=2.
%
% Return value:  -------------------------
%  linFWHH -- Full-width at half-height in linear units. 
%
% assert(size(linFWHH)==size(center))
% 
% Consider an arbitrary unimodal 1D tuning curve in linear scale, usually but
% not necessarily Gaussian. Introduce the following notation:
%   L: the left  point at half height
%   R: the right point at half height
%   C: the center of the Gaussian curve. The max height occurs for this point.
%   W: full-width at half-height in log scale (bandwidth)
%
% L and R satisfy the following simultaneous equations:
%  L+R = 2*C      % C is halfway between L and R
%  R/L = base^W   % the ratio between R and L is base^W
% Then:
%  L = 2*C/(base^W+1)
%  R = (base^W)*L
% The full width at half height FWHH on a linear scale:
%  FWHH = R-L = (base^W-1)*L = 2*C*(base^W-1)/(base^W+1)
%
% Compare with Equation 2.28 in Dayan & Abbott's (2001) book.
% Watson & Eckert (1994), JOSA, 11, 496-505, Eq. 2 give a similar formula.
%
% See also NormalDistrib_FWHH_to_StdDev, log2

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2013-06-18 TS: Wrote it


%% Handle input args
if (nargin<3) ; base = 2.0 ; end   % log is log2, so LogFWHH is in octaves


%% Do the real work
k = base^logFWHH;
linFWHH = 2 .* center .* (k-1) ./ (k+1);


%%% Return linFWHH
end %%% of file
