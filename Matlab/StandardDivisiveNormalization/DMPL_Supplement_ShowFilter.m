function [ output_args ] = DMPL_Supplement_ShowFilter(M,i_ori, i_spf)
%Shows image filters of the DNM
%
% Two image filters (Cos & Sin) specified by i_ori and i_spf are shown
%  Orientation (deg) = M.EarlyVis_spec.domain_orient_deg(i_ori);
%  Spatial-frequency (cpd) = M.EarlyVis_spec.domain_spFreq_cpd(i_spf);
%  Spatial-frequency (oct) = M.EarlyVis_spec.domain_spFreq_l2cpd(i_spf);

i_spf2 = FindClosestIdx(M.EarlyVis_spec.domain_channel_spFreq_cpd,...
                        M.EarlyVis_spec.domain_spFreq_cpd(i_spf));

i_maxX = size(M.stim_spec.gridX_deg,2);
i_maxY = size(M.stim_spec.gridY_deg,1);
rangeX = [M.stim_spec.gridX_deg(1,1), M.stim_spec.gridX_deg(1,i_maxX)];
rangeY = [M.stim_spec.gridY_deg(1,1), M.stim_spec.gridY_deg(i_maxY,1)];

maxVal = max(max(M.EarlyVis_filters.filters_Surr_CosImg{i_ori,i_spf2}(:)),...
             max(M.EarlyVis_filters.filters_Surr_SinImg{i_ori,i_spf2}(:)));
minVal = max(min(M.EarlyVis_filters.filters_Surr_CosImg{i_ori,i_spf2}(:)),...
             min(M.EarlyVis_filters.filters_Surr_SinImg{i_ori,i_spf2}(:)));

figure('Name','Filters');

subplot(1,2,1,'align');
imagesc(rangeX,rangeY,M.EarlyVis_filters.filters_Surr_CosImg{i_ori,i_spf2});
title('Cos filter');
%set(gca,'xtick',[],'ytick',[]);
axis(M.stim_spec.coordinatesStyle);
axis('image');
%truesize;
%set(gca,'visible','off');
caxis([minVal,maxVal]);
colorbar;

subplot(1,2,2,'align');
imagesc(rangeX,rangeY,M.EarlyVis_filters.filters_Surr_SinImg{i_ori,i_spf2});
title('Sin filter');
%set(gca,'xtick',[],'ytick',[]);
axis(M.stim_spec.coordinatesStyle);
axis('image');
%truesize;
%set(gca,'visible','off');
caxis([minVal,maxVal]);
colorbar;

output_args = sprintf(...
    'Orientation = %2.2f deg;  Spatial-frequency = %2.2f cpd (= %2.2f oct);',...
        M.EarlyVis_spec.domain_orient_deg(i_ori),...
        M.EarlyVis_spec.domain_channel_spFreq_cpd(i_spf2),...
        M.EarlyVis_spec.domain_channel_spFreq_l2cpd(i_spf2));

end

