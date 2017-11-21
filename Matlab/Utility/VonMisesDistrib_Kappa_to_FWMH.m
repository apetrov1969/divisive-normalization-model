function retval_fwmh_deg = VonMisesDistrib_Kappa_to_FWMH(kappa)

% VonMisesDistrib_Kappa_to_FWMH -- Convert the concentration parameter of von Mises distribution
%
% Conversion from Kappa (kurtosis) to FWMH (Full_Width at Middle Height) of von Mises distribution
% von Mises distribution: Exp(Kappa * cos(theta-Myu))/Constant, where Constant = 2*pi*Besseli(0,Kappa)
%
% Kappa: Concentration parameter of the von Mises distribution. Can be a matrix.
% FWMH (deg): Full Width at Middle Height (middle between min and max heights) of
%    von Mises distribution (-180~+180). Matrix of the same size as the input Kappa.
%
% FWMH [rad] = 2*acos(log(cosh(Kappa))/Kappa)
%
% For large Kappa (e.g., >100), the following approximation is used:
%   cosh(Kappa) = (exp(Kappa)+exp(-Kappa))/2  (approx)= Kappa - log(2)
%   FWMH [rad] = 2*acos(1 - log(2)/Kappa)
%
% For small Kappa (e.g., <.00001), Taylor-series approximation is used
%   log(cosh(Kappa))/Kappa = Kappa/2 - Kappa^3/12 + o(Kappa^5)
%
% Example 1:
%  kappa = 50;   % Concentration parameter of von Mises distribution
%  myu_deg = 30; % Center of von Mises distribution
%  fwmh_deg = VonMisesDistrib_Kappa_to_FWMH(kappa) % =19.1029 (deg)
%  hwmh_deg = fwmh_deg./2; % Half_Width at Middle_Height
%  height_at_peak = VonMisesDistribution(myu_deg,kappa,myu_deg)          % =2.8138
%  height_at_hwmh = VonMisesDistribution(myu_deg,kappa,myu_deg+hwmh_deg) % =1.4069
%  height_at_hwmh = VonMisesDistribution(myu_deg,kappa,myu_deg-hwmh_deg) % =1.4069
%  % height_at_peak is approximately equal to height_at_fwmh*2
%
% Example 2:
%  VonMisesDistrib_Kappa_to_FWMH([0 0.001 0.01 0.1 1 10 100 1000 10000 Inf])
%  ans = [180 179.9427 179.4270 174.2776 128.5845 42.9162 13.4999 4.2668 1.349 0]
%
% References:
%    Fisher, N. I. (1993) Statistical analysis of circular data. Cambridge University Press.
%
% See also VonMisesDistrib_FWMH_to_Kappa, VonMisesDistribuion, NormalDistrib_FWHH_to_StdDev.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.1 2013-06-20 AP: Added approximations for kappa<.0001 and kappa>100
% 1.0.0 2013-06-20 TS: Wrote it, from scratch


%% Calculate cos(hwmh)
cos_hwmh = NaN(size(kappa)) ;

% For kappa near zero, use the Taylor-series approximation
idx = find(kappa < .0001) ;
if (~isempty(idx))
    k = kappa(idx) ;
    cos_hwmh(idx) = k./2 - (k.^3)./12 ;
end

% For intermediate kappa values, use the exact equation
idx = find((kappa >= .0001) & (kappa < 100)) ;  
if (~isempty(idx))
    k = kappa(idx) ;
    cos_hwmh(idx) = log(cosh(k)) ./ k ;
end

% For large kappa, log(cosh(kappa)) is approximately equal to (kappa-log(2))
idx = find(kappa >= 100) ;  
if (~isempty(idx))
    cos_hwmh(idx) = 1 - log(2)./kappa(idx) ;
end


%% Calculate FWMH and return
retval_fwmh_deg = 2.*acos(cos_hwmh) .* (180/pi) ;

end
