function retval_kappa = VonMisesDistrib_FWMH_to_Kappa(fwmh_deg)

% VonMisesDistrib_FWMH_to_Kappa -- Convert the concentration parameter of von Mises distribution
%
% Conversion from FWMH (Full_Width at Middle Height) to Kappa (concentration) of von Mises distribution
% von Mises distribution: Exp(Kappa * cos(theta-mu))/Constant, where Constant = 2*pi*Besseli(0,Kappa)
%
% FWMH_deg -- Full Width at Middle Height (middele between min and max heights) of 
%   von Mises distribution (-180~+180). Must be a scalar. Matrices are not supported.
% Kappa -- Concentration parameter of the von Mises distribution. See Fisher (1993).
%
% Note that no deterministic transformation from FWMH to Kappa exists:
%   cos(FWMH/2) = log(cosh(Kappa))/Kappa          or equivalently
%   FWMH [rad] = 2*acos(log(cosh(Kappa))/Kappa)
%
% When cos(FWMH/2) is small, Kappa is large and cosh(Kappa) is approx = exp(Kappa)/2.
% Thus we have approximately:
%   cos(FWMH/2) = (Kappa - log(2))/Kappa = 1 - log(2)/Kappa
% Solving for Kappa, we get the approximation:
%   Kappa = log(2)/(1-cos(FWMH/2)), which is accurate to .0001 for FWMH<=55 deg.
%
% Example:
%    fwmh_deg = 22.3;       % Full_Width at Middle_Height
%    hwmh_deg = fwmh_deg/2; % Half_Width at Middle_Height
%    mu_deg = 45;
%    kappa = VonMisesDistrib_FWMH_to_Kappa(fwmh_deg) % =36.7670
%    height_at_peak = VonMisesDistribution(mu_deg,kappa,mu_deg)          % =2.4107
%    height_at_hwmh = VonMisesDistribution(mu_deg,kappa,mu_deg+hwmh_deg) % =1.2043
%    height_at_hwmh = VonMisesDistribution(mu_deg,kappa,mu_deg-hwmh_deg) % =1.2043
%    % height_at_peak is approximately equal to height_at_fwmh*2
%
% References:
%    Fisher, N. I. (1993) Statistical analysis of circular data. Cambridge University Press.
%
% See also VonMisesDistrib_Kappa_to_FWMH, VonMisesDistribuion, NormalDistrib_FWHH_to_StdDev.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.0 2013-06-20 TS & AP: Wrote it, from scratch


%% Enforce the range of FWMH
if (fwmh_deg<0)
    error('Full_Width at Middle_Height must be >=0');
    
elseif (fwmh_deg>180)
    warning('Full_Width at Middle_Height should be <=180');
    fprintf(1,'Kappa is set to be 0 (Uniform distribution)\n');
    retval_kappa = 0 ;
    return
end


%% Approximate formula -- accurate to .0001 for all FWMH<=55 deg.
if (fwmh_deg<=55)
    retval_kappa = log(2) / (1 - cos(fwmh_deg*pi/360)) ;
    return
end


%% Else FWMH > 55 deg and we have a transcendental equation
% We use a table computed in advance for every integer FWMH.
list_kappa = [...
      Inf,      NaN,      NaN,      NaN,      NaN, 505.7751, 371.6205, 284.5490, 224.8531, 182.1530,... %   1~ 10
 150.5597, 126.5305, 107.8301,  92.9919,  81.0212,  71.2240,  63.1044,  56.3001,  50.5416,  45.6251,... %  11~ 20
  41.3940,  37.7268,  34.5274,  31.7195,  29.2418,  27.0444,  25.0867,  23.3349,  21.7613,  20.3423,... %  21~ 30
  19.0584,  17.8931,  16.8320,  15.8632,  14.9763,  14.1622,  13.4132,  12.7226,  12.0845,  11.4936,... %  31~ 40
  10.9454,  10.4359,   9.9615,   9.5191,   9.1059,   8.7194,   8.3572,   8.0175,   7.6983,   7.3981,... %  41~ 50
   7.1154,   6.8489,   6.5973,   6.3595,   6.1346,   5.9216,   5.7198,   5.5282,   5.3464,   5.1735,... %  51~ 60
   5.0090,   4.8524,   4.7032,   4.5609,   4.4251,   4.2953,   4.1713,   4.0526,   3.9390,   3.8302,... %  61~ 70
   3.7258,   3.6257,   3.5295,   3.4371,   3.3483,   3.2628,   3.1805,   3.1012,   3.0248,   2.9511,... %  71~ 80
   2.8799,   2.8111,   2.7446,   2.6803,   2.6181,   2.5579,   2.4995,   2.4429,   2.3880,   2.3347,... %  81~ 90
   2.2829,   2.2326,   2.1837,   2.1361,   2.0898,   2.0447,   2.0008,   1.9580,   1.9162,   1.8754,... %  91~100
   1.8356,   1.7967,   1.7588,   1.7216,   1.6853,   1.6498,   1.6150,   1.5809,   1.5475,   1.5148,... % 101~110
   1.4827,   1.4512,   1.4203,   1.3900,   1.3602,   1.3310,   1.3022,   1.2739,   1.2461,   1.2188,... % 111~120
   1.1918,   1.1653,   1.1392,   1.1134,   1.0881,   1.0631,   1.0384,   1.0141,   0.9901,   0.9664,... % 121~130
   0.9429,   0.9198,   0.8970,   0.8744,   0.8521,   0.8300,   0.8082,   0.7866,   0.7652,   0.7440,... % 131~140
   0.7230,   0.7022,   0.6817,   0.6613,   0.6410,   0.6210,   0.6011,   0.5813,   0.5618,   0.5423,... % 141~150
   0.5230,   0.5038,   0.4848,   0.4658,   0.4470,   0.4283,   0.4097,   0.3912,   0.3728,   0.3545,... % 151~160
   0.3362,   0.3181,   0.3000,   0.2820,   0.2641,   0.2462,   0.2284,   0.2106,   0.1929,   0.1752,... % 161~170
   0.1576,   0.1400,   0.1224,   0.1049,   0.0873,   0.0699,   0.0524,   0.0349,   0.0175,   0];        % 171~180

fwmhRoundUp = ceil(fwmh_deg);
fwmhRoundDw = floor(fwmh_deg);
retval_kappa = list_kappa(fwmhRoundUp).*(fwmh_deg-fwmhRoundDw)...
       +list_kappa(fwmhRoundDw).*(fwmhRoundUp-fwmh_deg);

if fwmhRoundUp==fwmhRoundDw
    retval_kappa = list_kappa(fwmhRoundUp);
end

end
