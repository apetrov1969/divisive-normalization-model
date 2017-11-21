function full_width_half_height = NormalDistrib_StdDev_to_FWHH(stdev)

% NormalDistrib_StdDev_to_FWHH - Convert the bandwidth parameter
%
% full_width_half_height = NormalDistrib_StdDev_to_FWHH(stdev)
%
%FWHH_to_StdDev Conversion from FWHH to Standard Deviation (Sigma)
%   FWHH: Full Width at Half Height
%   Standard Deviation = sqrt(Variance) =Sigma
%   See Table 2.1 in Norma Graham's "Visual Pattern Analyzers" (1989, p. 49).
%   StDev = FWHH / 2sqrt(2ln2) = 0.425 FWHH;
%   FWHH  = StDev * 2sqrt(2ln2) = 2.35 StDev
%
% See also NormalDistrib_FWHH_to_StdDev, VonMisesDistrib_FWMH_to_Kappa

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2013-06-13 TS: Wrote it

full_width_half_height = stdev.*(2*sqrt(2*log(2)));

end

