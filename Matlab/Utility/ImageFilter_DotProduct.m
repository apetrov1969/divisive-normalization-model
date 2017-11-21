function  retval_responses = ImageFilter_DotProduct(filter,images)

%ImageFilter_DotProduct - Apply image filters using dot-product
%
% Inputs
%  filter: [y,x] or [y,x,filters]
%  images: [y,x] or [y,x,images]
%
% Outputs
%  retval_responses: [filters,images]
%

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2014-01-21 TS: Wrote it from scratch

size_filter = size(filter);
size_images = size(images);
dim_filter = size_filter(2);
dim_images = size_images(2);

if(size_filter(1)~=size_images(1) || size_filter(2)~=size_images(2))
    error('LinearModel: dimensions of inputs do not match.');
end

if(dim_images==2)
    num_images = 1;
else
    num_images = size_images(3);
end

if(dim_filter==2)
    temp_filter = reshape(filter,[1,size_filter(1)*size_filter(2)]); % [1,yx]
    temp_images = reshape(images,[size_images(1)*size_images(2),num_images]); % [yx,img]
    retval_responses = temp_filter * temp_images; % [1,img]
else
    num_filter = size_filter(3);
    temp_filter = reshape(filter,[size_filter(1)*size_filter(2),num_filter]); % [yx,filter]
    temp_images = reshape(images,[size_images(1)*size_images(2),num_images]); % [yx,img]
    retval_responses = temp_filter'*temp_images; % [filter,img]
end

%% Return M
end  %%% of file
