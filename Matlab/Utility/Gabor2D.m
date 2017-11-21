function retval_Gabor = Gabor2D(specs_gabor,gridX,gridY)

% Gabor2D -- Make a 2D Gabor image (sinusoidal grating under a Gaussian envelope)
%
% G = Gabor2D(specs_gabor,gridX,gridY)
%
% gridX, gridY: matrices indicating xy coordinates. The size of these
%               matrices must be identical to one another. The size of
%               Gabor image become the same with them.
%               Usually, gridX and gridY are produced by meshgrid for a
%               rectangular grid (see Example below).
%
% x = + j
% y = - i
% Clockwise rotation in the 'ij' coordinates.
% Counter-clockwise rotation in the 'xy' coordinates.
%               
%
% specs_gabor is a structure with fiels:
%   type:                 either 'Cos' or 'Sin'
%   amplitude:            grating is between +-amplitude
%   orientation_rad/_deg: orientatin of Gabor image
%                         'xy': counter-clock-wise
%                               x:left to right, index = 2
%                               y:bottom to top, index = 1
%                               0='||', 45='\\', 90='=', 135='//'
%                         'ij': clock-wise
%                               i:top to bottom, index = 1
%                               j:left to right, index = 2
%                               0='||', 45='//', 90='=', 135='\\'
%   frequency:            number of cycles of grating per unit
%   phase_rad/_deg:       phase of grating
%   centerX/Y:            center of grating and Gaussian envelope
%   peakHeight:           height of Gaussian envelope at peak (scalar value or 'normal')
%                          sum(Gauss(:))==1 if 'normal'
%   sigma_perpendicular:  sigma (width) of Gaussian envelope perpendicualr to grating
%   sigma_parallel:       sigma (width) of Gaussian envelope parallel to grating
%
% Example:
%    arrayX = [-100:+100];
%    arrayY = [-120:+120];
%    [gridX,gridY] = meshgrid(arrayX,arrayY); % rectangular grid
%    specs_gabor.type = 'Cos';
%    specs_gabor.amplitude = 0.5;
%    specs_gabor.orientation_deg = 10;
%    specs_gabor.frequency = 0.1;
%    specs_gabor.phase_rad = pi/4;
%    specs_gabor.centerX = 20;
%    specs_gabor.centerY = 60;
%    specs_gabor.peakHeight = 'normal';
%    specs_gabor.sigma_width  = 10;
%    specs_gabor.sigma_length = 30;
%    imgGabor = Gabor2D(specs_gabor,gridX,gridY) ;
%    imagesc(arrayX,arrayY,imgGabor); axis('xy')
%
% See also Gauss2D, Grating2D, meshgrid, gabor.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.2 2013-06-24 TS -- 'xy': Counter-clockwise, 'ij': Clockwise
% 1.0.1 2013-06-24 TS -- The coordinate style is consistently 'xy' with Clock-wise orientation starting from vertical
% 1.0.0 2013-06-17 TS -- Wrote it, based on utils/gabor.m


if ~isstruct(specs_gabor)
    error('Gabor2D.m: specs_gabor (1st input variable) must be a structure');
end
% if ~ismatrix(gridX)
%     error('Gabor2D.m: gridX (2nd input variable) must be a matrix');
% end
% if ~ismatrix(gridY)
%     error('Gabor2D.m: gridY (3rd input variable) must be a matrix');
% end
if ~isequal(size(gridX),size(gridY))
    error('Gabor2D.m: size of gridX and gridY must be identical to one another');
end

x=2;
y=1;


%% Set common specs for Gauss envelope and Grating pattern from specs_gabor
% Orientation of Gaussian envelope
if isfield(specs_gabor,'orientation_rad') && isfield(specs_gabor,'orientation_deg')
    error('Gabor2D.m: "orientation_rad" and "orientation_deg" are mutually exclusive to one another');
elseif isfield(specs_gabor,'orientation_rad')
    specs_gauss.orientation_rad = specs_gabor.orientation_rad;
elseif isfield(specs_gabor,'orientation_deg')
    specs_gauss.orientation_rad = specs_gabor.orientation_deg.*pi./180;
else
    specs_gauss.orientation_rad = 0;
end

% Center of Gaussian envelope
if isfield(specs_gabor,'center')
    if numel(specs_gabor.center)>=2
        specs_gauss.centerX = specs_gabor.center(x);
        specs_gauss.centerY = specs_gabor.center(y);
    else
        error('Gabor2D.m: center of Gabor patch should have two coordinates(y,x)');
    end
    if isfield(specs_gabor,'centerX') || isfield(specs_gabor,'centerY')
        error('Gabor2D.m: "center" and "centerX"/"centerY" are mutually exclusive to one another');
    end
elseif isfield(specs_gabor,'centerX') && isfield(specs_gabor,'centerY')
    specs_gauss.centerX = specs_gabor.centerX;
    specs_gauss.centerY = specs_gabor.centerY;
elseif isfield(specs_gabor,'centerX') || isfield(specs_gabor,'centerY')
    error('Gabor2D.m: center of Gabor patch should have both "centerX" and "centerY"');
else
    specs_gauss.centerX = 0;
    specs_gauss.centerY = 0;
end

% Copy from specs_gauss to specs_grating
specs_grating.orientation_rad = specs_gauss.orientation_rad;
specs_grating.centerX = specs_gauss.centerX;
specs_grating.centerY = specs_gauss.centerY;


%% Set specs for Gaussian envelope
% Height of Gaussian peak
if isfield(specs_gabor,'peakHeight') && isfield(specs_gabor,'gauss_peakHeight')
    error('Gabor2D.m: "gauss_peakHeight" and "peakHeight" are mutually exclusive to one another');
elseif isfield(specs_gabor,'peakHeight')
    specs_gauss.peakHeight = specs_gabor.peakHeight;
elseif isfield(specs_gabor,'gauss_peakHeight')
    specs_gauss.peakHeight = specs_gabor.gauss_peakHeight;
end

% Sigma (width) of Gaussian envelope
if isfield(specs_gabor,'sigma_width') && isfield(specs_gabor,'sigma_length')
    specs_gauss.sigmaX = specs_gabor.sigma_width;  % sigma along grating modulation:     (|) -> (||) -> (||...||)
    specs_gauss.sigmaY = specs_gabor.sigma_length; % sigma orthogonal to the modulation: (:) -> (==) -> (=.....=)
else
    error('Gabor2D.m: sigma_width and sigma_length are required');
end


%% Set specs for Grating pattern
% Type of grating
if isfield(specs_gabor,'type') && isfield(specs_gabor,'grating_type')
    error('Gabor2D.m: "type" and "grating_type" are mutually exclusive to one another');
elseif isfield(specs_gabor,'type')
    grating_type_str = specs_gabor.type;
elseif isfield(specs_gabor,'grating_type')
    grating_type_str = specs_gabor.grating_type;
else
    error('Gabor2D.m: grating_type (Cos/Sin) for Gabor patch is required');
end
str_cos={'cos','Cos','COS'};
str_sin={'sin','Sin','SIN'};
if nonzeros(strcmp(grating_type_str,str_cos))
    specs_grating.type='cos';
elseif nonzeros(strcmp(grating_type_str,str_sin))
    specs_grating.type='sin';
else
    error('Gabor2D.m: unknown grating_type (Cos/Sin) for Gabor patch');
end

% Amplitude of grating
if isfield(specs_gabor,'amplitude') && isfield(specs_gabor,'grating_amplitude')
    error('Gabor2D.m: "amplitude" and "grating_amplitude" are mutually exclusive to one another');
elseif isfield(specs_gabor,'amplitude')
    specs_grating.amplitude = specs_gabor.amplitude;
elseif isfield(specs_gabor,'grating_amplitude')
    specs_grating.amplitude = specs_gabor.grating_amplitude;
else
    specs_grating.amplitude = 1;
end

% Frequency of grating
if isfield(specs_gabor,'frequency') && isfield(specs_gabor,'grating_frequency')
    error('Gabor2D.m: "frequency" and "grating_frequency" are mutually exclusive to one another');
elseif isfield(specs_gabor,'frequency')
    specs_grating.frequency = specs_gabor.frequency;
elseif isfield(specs_gabor,'grating_frequency')
    specs_grating.frequency = specs_gabor.grating_frequency;
else
    error('Gabor2D.m: frequency for Gabor patch is required');
end

% Phase of grating
str_phase={'phase_rad','phase_deg','grating_phase_rad','grating_phase_deg'};
if sum(isfield(specs_gabor,str_phase))>1
    error('Gabor2D.m: "phase_rad/_deg" and "grating_phase_rad/_deg" are mutually exclusive to one another');
elseif isfield(specs_gabor,'phase_rad')
    specs_grating.phase_rad = specs_gabor.phase_rad;
elseif isfield(specs_gabor,'grating_phase_rad')
    specs_grating.phase_rad = specs_gabor.grating_phase_rad;
elseif isfield(specs_gabor,'phase_deg')
    specs_grating.phase_rad = specs_gabor.phase_deg.*pi./180;
elseif isfield(specs_gabor,'grating_phase_deg')
    specs_grating.phase_rad = specs_gabor.grating_phase_deg.*pi./180;
else
    specs_grating.phase_rad = 0;
end


%% Make Gabor image (sinusoidal grating under a Gaussian envelope)
gaussEnvelope  = Gauss2D(specs_gauss,gridX,gridY);
gratingPattern = Grating2D(specs_grating,gridX,gridY);

retval_Gabor = gratingPattern .* gaussEnvelope;

end


    