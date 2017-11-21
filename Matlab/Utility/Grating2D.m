function retval_Grating = Grating2D(specs_grating,gridX,gridY)

% Grating2D -- Make a sinusoidal grating (grayscale image)
%
% G = Grating2D(specs_grating,gridX,gridY)
%
% gridX, gridY: matrices indicating xy coordinates. The size of these
%               matrices must be identical to one another. The size of
%               Grating image become the same with them.
%               Usually, gridX and gridY are produced by meshgrid for a
%               rectangular grid (see Example below).
%
% specs_grating is a structure with fields
%   type:                 either 'Cos' or 'Sin'
%   amplitude:            grating is between +-amplitude
%   orientation_rad/_deg: orientatin of grating
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
%   centerX/Y:            center of grating
%
% Example:
%    arrayX = [-100:+100];
%    arrayY = [-120:+120];
%    [gridX,gridY] = meshgrid(arrayX,arrayY); % rectangular grid
%    specs_grating.type = 'Cos';
%    specs_grating.amplitude = 0.5;
%    specs_grating.orientation_rad = pi/18;
%    specs_grating.frequency = 0.05;
%    specs_grating.phase_deg = 45;
%    specs_grating.centerX = 0;
%    specs_grating.centerY = 0;
%    imgGrating = Grating2D(specs_grating,gridX,gridY) ;
%    imagesc(arrayX,arrayY,imgGrating); axis('xy')
%
% See also Gabor2D, Gauss2D, meshgrid, grating.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.3 2013-06-24 TS -- 'xy': Counter-clockwise, 'ij': Clockwise
% 1.0.2 2013-06-27 TS -- Positive rotation refers rotation from the x-axis to the y-axis in any coordinate setting
% 1.0.1 2013-06-24 TS -- The coordinate style is consistently 'xy' with Clock-wise orientation starting from vertical
% 1.0.0 2013-06-17 TS -- Wrote it, based on utils/grating.m


if ~isstruct(specs_grating)
    error('Grating2D.m: specs_grating (1st input variable) must be a structure');
end
% if ~ismatrix(gridX)
%     error('Grating2D.m: gridX (2nd input variable) must be a matrix');
% end
% if ~ismatrix(gridY)
%     error('Grating2D.m: gridY (3rd input variable) must be a matrix');
% end
if ~isequal(size(gridX),size(gridY))
    error('Grating2D.m: size of gridX and gridY must be identical to one another');
end

x=2;
y=1;


%% Get specs
if isfield(specs_grating,'type')
    type_str = specs_grating.type;
else
    error('Grating2D.m: type (Cos/Sin) is required');
end
str_cos={'cos','Cos','COS'};
str_sin={'sin','Sin','SIN'};
if nonzeros(strcmp(type_str,str_cos))
    type='cos';
elseif nonzeros(strcmp(type_str,str_sin))
    type='sin';
else
    error('Grating2D.m: unknown type (Cos/Sin)');
end

% Amplitude of grating
if isfield(specs_grating,'amplitude')
    amplitude = specs_grating.amplitude;
else
    amplitude = 1;
end

% Frequency of grating
if isfield(specs_grating,'frequency')
    frequency = specs_grating.frequency;
else
    error('Grating2D.m: frequency of grating is required');
end

% Phase of grating
if isfield(specs_grating,'phase_rad')
    phase_rad = specs_grating.phase_rad;
elseif isfield(specs_grating,'phase_deg')
    phase_rad = specs_grating.phase_deg.*pi./180;
else
    phase_rad = 0;
end

% Orientation of grating
if isfield(specs_grating,'orientation_rad') && isfield(specs_grating,'orientation_deg')
    error('Grating2D.m: "orientation_rad" and "orientation_deg" are mutually exclusive to one another');
elseif isfield(specs_grating,'orientation_rad')
    orientation_rad = specs_grating.orientation_rad;
elseif isfield(specs_grating,'orientation_deg')
    orientation_rad = specs_grating.orientation_deg.*pi./180;
else
    orientation_rad = 0;
end

% Center of grating
if isfield(specs_grating,'center')
    if numel(specs_grating.center)>=2
        centerX = specs_grating.center(x);
        centerY = specs_grating.center(y);
    else
        error('Grating2D.m: center of grating should have two coordinates(y,x)');
    end
    if isfield(specs_grating,'centerX') || isfield(specs_grating,'centerY')
        error('Grating2D.m: "center" and "centerX"/"centerY" are mutually exclusive to one another');
    end
elseif isfield(specs_grating,'centerX') && isfield(specs_grating,'centerY')
    centerX = specs_grating.centerX;
    centerY = specs_grating.centerY;
elseif isfield(specs_grating,'centerX') || isfield(specs_grating,'centerY')
    error('Grating2D.m: center of grating should have both "centerX" and "centerY"');
else
    centerX = 0;
    centerY = 0;
end


%% Rotate grid for -orientation (== rotate picture for +orientation)
sinOri = sin(-orientation_rad);
cosOri = cos(-orientation_rad);
rotX = (gridX-centerX).*cosOri - (gridY-centerY).*sinOri;
rotY = (gridX-centerX).*sinOri + (gridY-centerY).*cosOri;

if strcmp(type,'cos')
    retval_Grating = cos(rotX.*2.*pi.*frequency-phase_rad) .* amplitude;
elseif strcmp(type,'sin')
    retval_Grating = sin(rotX.*2.*pi.*frequency-phase_rad) .* amplitude;
end

end
