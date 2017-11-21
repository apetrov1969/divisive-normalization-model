function [ output_args ] = DMPL_Supplement_ShowLocationsInImage(M, image)
%Shows locations of receptive-field center superimposed on an image
%
% Locations of receptive-field centers of the DNM neurons are shown
% on a given image
%  output_args = locations (y,x);

if isempty(image); image = zeros(size(M.stim_spec.gridX_deg)); end;
if ismatrix(image)==0; error('Invalid input image (1)'); end;
if max(image(:))>1; error('Invalid input image (2)'); end;

locations = M.EarlyVis_spec.domain_rfLoc_deg;

i_maxX = size(image,2);
i_maxY = size(image,1);
rangeX = [M.stim_spec.gridX_deg(1,1), M.stim_spec.gridX_deg(1,i_maxX)];
rangeY = [M.stim_spec.gridY_deg(1,1), M.stim_spec.gridY_deg(i_maxY,1)];

figure('Name','RF centers');
imagesc(rangeX,rangeY,image(:,:));
%set(gca,'xtick',[],'ytick',[]);
axis(M.stim_spec.coordinatesStyle);
axis('image');
truesize;
%set(gca,'visible','off');
colormap('gray');
caxis([0,1]);
hold on;

num_locations = size(locations,1);
for l = 1:num_locations
    plot(locations(l,2),locations(l,1),'r.','MarkerSize',20);
    text(locations(l,2),locations(l,1),num2str(l),'Color','Green');
end

hold off;

output_args = locations;

end

