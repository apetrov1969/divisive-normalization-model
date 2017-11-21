function [retval_y,retval_x,retval_num] = FindLocalExtrema2D(image2D)

% FindLocalExtrema2D -- Make a sinusoidal grating (grayscale image)
%
% G = Grating2D(specs_grating,gridX,gridY)

image2D = squeeze(image2D);

if ~ismatrix(image2D)
    error('Grating2D.m: gridX (2nd input variable) must be a matrix');
end

y=1;
x=2;

size_img = size(image2D);

paddedImage = zeros(size_img(y)+2,size_img(x)+2);
regionY = (1:size_img(y))+1;
regionX = (1:size_img(x))+1;
paddedImage = paddedImage+min(image2D(:))-1;
paddedImage(regionY,regionX) = image2D;

retval_num = 0;
extrema=zeros(size_img);
for iy=1:size_img(y)
    for ix=1:size_img(x)
        if paddedImage(iy+1,ix+1) > paddedImage(iy+1-1,ix+1)...   % |
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1+1,ix+1)...   % |
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1,  ix+1-1)... % -
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1,  ix+1+1)... % -
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1-1,ix+1-1)... % /
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1+1,ix+1+1)... % /
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1+1,ix+1-1)... % \
        && paddedImage(iy+1,ix+1) > paddedImage(iy+1-1,ix+1+1)    % \
            retval_num = retval_num+1;
            extrema(iy,ix)=1;
        end
    end % for x=2:maxX-1
end % for y=2:maxY-1

[retval_y,retval_x] = find(extrema);

end