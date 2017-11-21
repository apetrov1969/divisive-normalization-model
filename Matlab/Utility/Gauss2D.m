function retval_Gauss = Gauss2D(specs_gauss,gridX,gridY)

% Gauss2D -- Make a 2D grayscale image of a Gaussian blob
%
% G = Gauss2D(specs_gauss,gridX,gridY)
%
% gridX, gridY: matrices indicating xy coordinates. The size of these
%               matrices must be identical to one another. The size of
%               the Gaussian image become the same with them.
%               Usually, gridX and gridY are produced by meshgrid for a
%               rectangular grid (see Example below).
%
% specs_gauss is a structure with fiels:
% orientation_rad/_deg: orientatin of Gabor image
%                         'xy': counter-clock-wise
%                               x:left to right, index = 2
%                               y:bottom to top, index = 1
%                               0='||', 45='\\', 90='=', 135='//'
%                         'ij': clock-wise
%                               i:top to bottom, index = 1
%                               j:left to right, index = 2
%                               0='||', 45='//', 90='=', 135='\\'
% centerX/Y:            center of Gaussian distribution
% peakHeight:           height of Gaussian distribution at peak (scalar value or 'normal')
% sigmaX:               sigma (width) of Gaussian distribution along the x-axis if orientation==0
% sigmaY:               sigma (width) of Gaussian distribution along the y-axis if orientation==0
%
% Example:
%    arrayX = [-100:+100];
%    arrayY = [-120:+120];
%    [gridX,gridY] = meshgrid(arrayX,arrayY); % rectangular grid
%    specs_gauss.orientation_deg = 10;
%    specs_gauss.centerX = 20;
%    specs_gauss.centerY = 60;
%    specs_gauss.peakHeight = 'normal';
%    specs_gauss.sigmaX = 10;
%    specs_gauss.sigmaY = 30;
%    imgGauss = Gauss2D(specs_gauss,gridX,gridY) ;
%    imagesc(arrayX,arrayY,imgGauss); axis('xy')
%
% See also Gabor2D, Grating2D, meshgrid, normpdf.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu

% 1.0.2 2013-06-24 TS -- 'xy': Counter-clockwise, 'ij': Clockwise
% 1.0.1 2013-06-24 TS -- The coordinate style is consistently 'xy' with Clock-wise orientation starting from vertical
% 1.0.0 2013-06-17 TS -- Wrote it, based on utils/gabor.m


if ~isstruct(specs_gauss)
    error('Gauss2D.m: specs_gauss (1st input variable) must be a structure');
end
% if ~ismatrix(gridX)
%     error('Gauss2D.m: gridX (2nd input variable) must be a matrix');
% end
% if ~ismatrix(gridY)
%     error('Gauss2D.m: gridY (3rd input variable) must be a matrix');
% end
if ~isequal(size(gridX),size(gridY))
    error('Gauss2D.m: size of gridX and gridY must be identical to one another');
end

x=2;
y=1;


%% Get specs
% Height of Gaussian distribution at peak
if isfield(specs_gauss,'peakHeight')
    if ischar(specs_gauss.peakHeight)
        if strcmp(specs_gauss.peakHeight,'normal') || strcmp(specs_gauss.peakHeight,'Normal') || strcmp(specs_gauss.peakHeight,'NORMAL')
            peakHeight = 0;
            normalize = true;
        else
            error('Gauss2D.m: peakHeight of Gaussian distribution is either scalar or "normal"');
        end
    else
        peakHeight = specs_gauss.peakHeight;
        normalize = false;
    end
else
    peakHeight = 0;
    normalize = true;
end

% Orientation of Gaussian distribution
if isfield(specs_gauss,'orientation_rad') && isfield(specs_gauss,'orientation_deg')
    error('Gauss2D.m: "orientation_rad" and "orientation_deg" are mutually exclusive to one another');
elseif isfield(specs_gauss,'orientation_rad')
    orientation_rad = specs_gauss.orientation_rad;
elseif isfield(specs_gauss,'orientation_deg')
    orientation_rad = specs_gauss.orientation_deg.*pi./180;
else
    orientation_rad = 0;
end

% Sigma (width) of Gaussian distribution
if isfield(specs_gauss,'sigma')
    if numel(specs_gauss.sigma)==1
        sigmaX = specs_gauss.sigma;
        sigmaY = specs_gauss.sigma;
    else
        sigmaX = specs_gauss.sigma(x);
        sigmaY = specs_gauss.sigma(y);
    end
    if isfield(specs_gauss,'sigmaX') || isfield(specs_gauss,'sigmaY')
        error('Grating2D.m: "sigma" and "sigmaX"/"sigmaY" are mutually exclusive to one another');
    end
elseif isfield(specs_gauss,'sigmaX') && isfield(specs_gauss,'sigmaY')
    sigmaX = specs_gauss.sigmaX; % stdev along the x-axis if orientation==0
    sigmaY = specs_gauss.sigmaY; % stdev along the y-axis if orientation==0
elseif isfield(specs_gauss,'sigmaX') || isfield(specs_gauss,'sigmaY')
    error('Grating2D.m: sigma of Gaussian distribution shoul have both "sigmaX" and "sigmaY" ');
else
    error('Gauss2D.m: sigma of Gaussian distribution is required');
end

% Center of Gaussian distribution
if isfield(specs_gauss,'center')
    if numel(specs_gauss.center)>=2
        centerX = specs_gauss.center(x);
        centerY = specs_gauss.center(y);
    else
        error('Gauss2D.m: center of Gaussian distribution should have two coordinates(y,x)');
    end
    if isfield(specs_gauss,'centerX') || isfield(specs_gauss,'centerY')
        error('Gauss2D.m: "center" and "centerX"/"centerY" are mutually exclusive to one another');
    end
elseif isfield(specs_gauss,'centerX') && isfield(specs_gauss,'centerY')
    centerX = specs_gauss.centerX;
    centerY = specs_gauss.centerY;
elseif isfield(specs_gauss,'centerX') || isfield(specs_gauss,'centerY')
    error('Gauss2D.m: center of grating should have both "centerX" and "centerY"');
else
    centerX = 0;
    centerY = 0;
end


%% Rotate grid for -orientation (== rotate picture for +orientation))
sinOri = sin(-orientation_rad);
cosOri = cos(-orientation_rad);
rotX = (gridX-centerX).*cosOri - (gridY-centerY).*sinOri;
rotY = (gridX-centerX).*sinOri + (gridY-centerY).*cosOri;


%% Gaussian envelope (rotX and rotY are independent)
gaussX = exp(-0.5 .* (rotX./sigmaX).^2);
gaussY = exp(-0.5 .* (rotY./sigmaY).^2);
gaussXY = gaussX .* gaussY;

if normalize
    retval_Gauss = gaussXY / sum(gaussXY(:));
else
    retval_Gauss = gaussXY * peakHeight;
end

end