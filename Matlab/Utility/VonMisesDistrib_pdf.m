function retval_vonMises = VonMisesDistrib_pdf(myu_deg,kappa,theta_deg)

% VonMisesDistrib_pdf -- Probability density function of the von Mises circular distrib
%
% The von Mises distribution is defined by the equation: 
%  Exp(Kappa * cos(theta-Myu))/Constant, where 
%  Constant = 2*pi*Besseli(0,Kappa)
%
% Inputs:   %%%%%%%%%%%%%%
%   myu_deg: Center of distribution (0~360 deg), scalar
%   kappa: Concentration parameter of the von Mises distribution, scalar
%   theta_deg: Grid of points over which the return values are calcualted.
%
% Return value:  %%%%%%%%%%
%   vonMises: probability density values. Matrix of the same size as theta_grid
%
% Example:
%    theta = 0:359; % deg
%    myu = 45; % deg
%    kappa = 10;
%    vMpdf = VonMisesDistrib_pdf(myu,kappa,theta);
%    plot(theta,vMpdf) ; sum(vMpdf)./360.*(2*pi) % =1
%
% References:
%    Fisher, N. I. (1993) Statistical analysis of circular data. Cambridge Univ Press.
%
% See also VonMisesDistrib_FWMH_to_Kappa, VonMisesDistrib_Kappa_to_FWMH, pdf.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.1 2013-06-21 AP: Return the delta-function pdf for Kappa>700
% 1.0.0 2013-06-20 TS: Wrote it, from scratch

if (kappa<=700)
    myu_rad   =   myu_deg .* (pi/180) ;
    theta_rad = theta_deg .* (pi/180) ;

    retval_vonMises = exp(kappa .* cos(theta_rad-myu_rad)) / (2*pi*besseli(0,kappa));

else  % Kappa tends to Inf
    % exp(710)       = inf; (MATLAB R2012b)
    % besseli(0,701) = inf; (MATLAB R2012b)
    
    % Return the delta-function pdf
    retval_vonMises = zeros(size(theta_deg)) ;
    retval_vonMises(theta_deg==myu_deg) = Inf ;
end


%%% Return retval_vomMises
end %%% of file
