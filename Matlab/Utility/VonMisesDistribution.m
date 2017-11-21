function retval_vonMises = VonMisesDistribution(myu_deg,kappa,theta_deg)
% von Mises distribution: Exp(Kappa * cos(theta-Myu))/Constant, where Constant = 2*pi*Besseli(0,Kappa)
% myu_deg: Center of distribution (0~360 deg)
% kappa: Kurtosis, characterizing a ratio between Min and Max of von Mises distribution
%
% Example:
%    theta = 1:360; % deg
%    myu = 45; % deg
%    kappa = 100;
%    distribution = VonMisesDistribution(myu,kappa,theta);
%    sum(distribution)./360.*(2*pi) % =1
%
% References:
%    Fisher, N. I. (1993) Statistical analysis of circular data. Cambridge University Press.
%
% See also VonMisesDistrib_FWMH_to_Kappa, VonMisesDistrib_Kappa_to_FWMH.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.0 2013-06-20 TS: Wrote it, from scratch

if(kappa>700)
    % exp(710)       = inf; (MATLAB R2012b)
    % besseli(0,701) = inf; (MATLAB R2012b)
    error('Kappa should be <=700 becase of overflow');
end

myu   =   myu_deg./180.*pi;
theta = theta_deg./180.*pi;
retval_vonMises = exp(kappa .* cos(theta-myu)) / (2*pi*besseli(0,kappa));

end
