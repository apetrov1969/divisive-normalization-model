function standard_deviation = NormalDistrib_FWHH_to_StdDev(fwhh)

% NormalDistrib_FWHH_to_StdDev - Convert the bandwidth parameter of the normal distrib
%
% sigma = NormalDistrib_FWHH_to_StdDev(fwhh)
%
%   FWHH:  Full Width at Half Height
%   Sigma: Standard Deviation = sqrt(Variance)
%
%   See Table 2.1 in Norma Graham's "Visual Pattern Analyzers" (1989, p. 49).
%   StDev = FWHH / 2sqrt(2ln2) = 0.425 FWHH;
%   FWHH  = StDev * 2sqrt(2ln2) = 2.35 StDev
%
% See also NormalDistrib_StdDev_to_FWHH, Log_to_Linear_FWHH, VonMisesDistrib_FWMH_to_Kappa

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2013-06-18 TS: Wrote it

standard_deviation=fwhh./(2*sqrt(2*log(2)));

end

