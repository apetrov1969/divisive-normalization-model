function retval_image2D = SquareWave2D(specs_wave,gridX,gridY)

% SquareWave2D -- Make a sinusoidal wave (grayscale image)
%
% retval_image2D = SquareWave2D(specs_wave,gridX,gridY)
%
% gridX, gridY: matrices indicating xy coordinates. The size of these
%               matrices must be identical to one another. The size of
%               Wave image become the same with them.
%               Usually, gridX and gridY are produced by meshgrid for a
%               rectangular grid (see Example below).
%
% specs_wave is a structure with fields
%   type:                 either 'Cos' or 'Sin'
%   amplitude:            wave is between +-amplitude
%   orientation_rad/_deg: orientatin of grating
%                         'xy': counter-clock-wise
%                               x:left to right, index = 2
%                               y:bottom to top, index = 1
%                               0='||', 45='\\', 90='=', 135='//'
%                         'ij': clock-wise
%                               i:top to bottom, index = 1
%                               j:left to right, index = 2
%                               0='||', 45='//', 90='=', 135='\\'
%   frequency:            number of cycles of wave per unit
%   phase_rad/_deg:       phase of wave
%   centerX/Y:            center of wave
%
% Example:
%    arrayX = [-100:+100];
%    arrayY = [-120:+120];
%    [gridX,gridY] = meshgrid(arrayX,arrayY); % rectangular grid
%    specs_wave.type = 'Cos';
%    specs_wave.amplitude = 0.5;
%    specs_wave.orientation_rad = pi/18;
%    specs_wave.frequency = 0.05;
%    specs_wave.phase_deg = 45;
%    specs_wave.centerX = 0;
%    specs_wave.centerY = 0;
%    imgWave = SquareWave2D(specs_wave,gridX,gridY) ;
%    imagesc(arrayX,arrayY,imgWave); axis('xy')
%
% See also Grating2D, meshgrid, wave.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.1 2013-06-24 TS -- 'xy': Counter-clockwise, 'ij': Clockwise
% 1.0.0 2014-01-23 TS -- Wrote it, based on Grating2D.m


if ~isstruct(specs_wave)
    error('SquareWave2D.m: specs_wave (1st input variable) must be a structure');
end
if ~ismatrix(gridX)
    error('SquareWave2D.m: gridX (2nd input variable) must be a matrix');
end
if ~ismatrix(gridY)
    error('SquareWave2D.m: gridY (3rd input variable) must be a matrix');
end
if ~isequal(size(gridX),size(gridY))
    error('SquareWave2D.m: size of gridX and gridY must be identical to one another');
end

x=2;
y=1;


%% Get specs
if isfield(specs_wave,'type')
    type_str = specs_wave.type;
else
    error('SquareWave2D.m: type (Cos/Sin) is required');
end
str_cos={'cos','Cos','COS'};
str_sin={'sin','Sin','SIN'};
if nonzeros(strcmp(type_str,str_cos))
    type='cos';
elseif nonzeros(strcmp(type_str,str_sin))
    type='sin';
else
    error('SquareWave2D.m: unknown type (Cos/Sin)');
end

% Amplitude of wave
if isfield(specs_wave,'amplitude')
    amplitude = specs_wave.amplitude;
else
    amplitude = 1;
end

% Frequency of wave
if isfield(specs_wave,'frequency')
    frequency = specs_wave.frequency;
else
    error('SquareWave2D.m: frequency of wave is required');
end

% Phase of wave
if isfield(specs_wave,'phase_rad')
    phase_rad = specs_wave.phase_rad;
elseif isfield(specs_wave,'phase_deg')
    phase_rad = specs_wave.phase_deg.*pi./180;
else
    phase_rad = 0;
end

% Orientation of wave
if isfield(specs_wave,'orientation_rad') && isfield(specs_wave,'orientation_deg')
    error('SquareWave2D.m: "orientation_rad" and "orientation_deg" are mutually exclusive to one another');
elseif isfield(specs_wave,'orientation_rad')
    orientation_rad = specs_wave.orientation_rad;
elseif isfield(specs_wave,'orientation_deg')
    orientation_rad = specs_wave.orientation_deg.*pi./180;
else
    orientation_rad = 0;
end

% Center of wave
if isfield(specs_wave,'center')
    if numel(specs_wave.center)>=2
        centerX = specs_wave.center(x);
        centerY = specs_wave.center(y);
    else
        error('SquareWave2D.m: center of wave should have two coordinates(y,x)');
    end
    if isfield(specs_wave,'centerX') || isfield(specs_wave,'centerY')
        error('SquareWave2D.m: "center" and "centerX"/"centerY" are mutually exclusive to one another');
    end
elseif isfield(specs_wave,'centerX') && isfield(specs_wave,'centerY')
    centerX = specs_wave.centerX;
    centerY = specs_wave.centerY;
elseif isfield(specs_wave,'centerX') || isfield(specs_wave,'centerY')
    error('SquareWave2D.m: center of wave should have both "centerX" and "centerY"');
else
    centerX = 0;
    centerY = 0;
end

% Type of wave (modifying phase_rad)
if strcmp(type,'cos')
    phase_rad = phase_rad-pi./2;
% elseif strcmp(type,'sin')
end


%% Rotate grid for -orientation (== rotate picture for +orientation)
sinOri = sin(-orientation_rad);
cosOri = cos(-orientation_rad);
rotX = (gridX-centerX).*cosOri - (gridY-centerY).*sinOri;
rotY = (gridX-centerX).*sinOri + (gridY-centerY).*cosOri;

retval_image2D = square(rotX.*2.*pi.*frequency-phase_rad) .* amplitude;

end
