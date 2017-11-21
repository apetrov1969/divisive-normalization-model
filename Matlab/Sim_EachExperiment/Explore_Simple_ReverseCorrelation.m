function [retval_mainResults,...
          retval_allResults] = Explore_Simple_ReverseCorrelation(M, specs_simulation)

%Explore_Simple_ReverseCorrelation - Measuring simple-cell receptive field using reverse correlation method
%
% Visual stimuli:
%  Random noise pattern
%
% References:
%  Gardner, Anzai, Ohzawa, & Freeman (1999)
%  Nishimoto, Ishida, & Ohzawa (2006)
%  Moore IV & Freeman (2012)
%  Ringach (2002)
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2); %#ok
maxGridR_deg = max(abs(gridX(:))); %#ok

num_ev_ori  = evS.num_orient; %#ok
num_ev_spf  = evS.num_spFreq; %#ok
num_ev_type = evS.num_cellType; % [1 Complex cell, 4 Simple cells]


%% Set up simulation
%- Get indices of neurons for analysing the results
neuron_orient_deg = specs_simulation.neuron_orient_deg;
neuron_spFreq_cpd = specs_simulation.neuron_spFreq_cpd;
idx_neuron_orient = FindClosestIdx(evS.domain_orient_deg,neuron_orient_deg);
idx_neuron_spFreq = FindClosestIdx(evS.domain_spFreq_cpd,neuron_spFreq_cpd);
idx_neuron_rfLoc = 1;

%- Set specs of stimuli
specs_images.size_block_pix = specs_simulation.size_block_pix;
specs_images.num_images = specs_simulation.num_images;

specs_images.contrast = 0.9; % Gardner, Anazai, Ohzawa & Freeman (1999)

numConditions = specs_images.num_images;

noise_block_pix   = specs_images.size_block_pix;
noise_image_block = floor(stS.imageSize_pix./noise_block_pix);
setStimuli = NaN([stS.imageSize_pix, numConditions]);
for cf=1:numConditions
    tempImg = rand(noise_image_block);
    tempImg = tempImg*stS.rangeLum - stS.bgLum; % -rangeLum/2~+rangeLum/2
    tempImg = tempImg*specs_images.contrast + stS.bgLum;
    setStimuli(:,:,cf) = imresize(tempImg,[stS.imageSize_pix],'nearest');
end


%% Apply Dimple Early Vision process
[resp_CentComp,...
 resp_CentSimp,...
 resp_SuppChan] = DMPL_EarlyVis_PoolRectifiedLinChannels(M, setStimuli);
[resp_DivNorm,...
 resp_SuppDrv,...
 resp_StimDrv]  = DMPL_EarlyVis_MeanChannelResp(M,resp_CentComp,...
                                                  resp_CentSimp,...
                                                  resp_SuppChan);          

%- Set 'retval_allResults'
if(nargout>=2)
    retval_allResults.resp_CentComp = resp_CentComp;
    retval_allResults.resp_CentSimp = resp_CentSimp;
    retval_allResults.resp_SuppChan = resp_SuppChan;
    retval_allResults.resp_DivNorm  = resp_DivNorm;
    retval_allResults.resp_SuppDrv  = resp_SuppDrv;
    retval_allResults.resp_StimDrv  = resp_StimDrv;
end


%% Analyze results
%- Organize results
results = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,img]
clear resp_CentComp resp_CentSimp resp_SuppChan;

%- Derive receptive fields
receptiveFields = zeros([stS.imageSize_pix,num_ev_type]);
for tp=1:num_ev_type
    for im=1:numConditions
        receptiveFields(:,:,tp) = receptiveFields(:,:,tp) + results(tp,im).*setStimuli(:,:,im)/1000;
    end
end

%- Analyze derived receptive fields in Fourier domain
background = squeeze(mean(mean(receptiveFields,1),2));
for tp=1:num_ev_type
    receptiveFields(:,:,tp) = receptiveFields(:,:,tp)-background(tp);
end

size_imgPadded = M.EarlyVis_filters.size_imgPadded;
size_padRegion = M.EarlyVis_filters.size_leftBtmPad;
regionY = (1:stS.imageSize_pix(1))+size_padRegion(1); %y: 1st index
regionX = (1:stS.imageSize_pix(2))+size_padRegion(2); %x: 2nd index
paddedRF  = zeros(size_imgPadded);
fourierRF = NaN([size_imgPadded,num_ev_type]);
for tp=1:num_ev_type
    paddedRF(:) = 0;
    paddedRF(regionY,regionX) = receptiveFields(:,:,tp);
    fourierRF(:,:,tp) = fftshift(abs(fft2(paddedRF)));
end


%% Derive Gauss distributions from Fourier receptive fields
width_imgPadded_deg = size_imgPadded(1)*stS.degPerPixel; % Assume imgPadded is square
oyFourier = size_imgPadded(1)/2+1;
oxFourier = size_imgPadded(2)/2+1;
domainX = (1:size_imgPadded(2));
domainY = (1:size_imgPadded(1));
[gridX,gridY] = meshgrid(domainX,domainY); % deg

gauss1D = @(grid,sig,cen) ( exp(-(grid-cen).^2 ./ (2*sig.^2)) );

pairGauss2D = @(sigY,sigX,cenX) ( gauss1D(gridY,sigY/10,oyFourier) .* gauss1D(gridX,sigX/10,oxFourier+cenX)...
                                 +gauss1D(gridY,sigY/10,oyFourier) .* gauss1D(gridX,sigX/10,oxFourier-cenX));

fittedFourierRF = NaN([size_imgPadded,num_ev_type]);
opt_values=cell(5,1);

temp_cx = neuron_spFreq_cpd*width_imgPadded_deg;
temp_sx = 10*width_imgPadded_deg*M.EarlyVis_filters.sigma_clrf_spFreq_cpd(idx_neuron_spFreq);
temp_sy = 10*width_imgPadded_deg*M.EarlyVis_filters.sigma_clrf_orient_rad * evS.domain_spFreq_cpd(idx_neuron_spFreq);
temp_values = [temp_sy,temp_sx,temp_cx,1];

for tp=2:num_ev_type
    temp_values(4) = Max2D(squeeze(fourierRF(:,:,tp)));
    
    sum_sqdiff = @(v) mean(mean((v(4)*pairGauss2D(v(1),v(2),v(3)) - fourierRF(:,:,tp)).^2,1),2);
    opt_values{tp} = fminsearch(sum_sqdiff,temp_values);
    
    fittedFourierRF(:,:,tp) = opt_values{tp}(4) * pairGauss2D(opt_values{tp}(1),opt_values{tp}(2),opt_values{tp}(3));
    
    opt_values{tp}(1) = opt_values{tp}(1)/10;
    opt_values{tp}(2) = opt_values{tp}(2)/10;
end


%% Derive bandwidths of orientation/spatial-frequency
derived_spf_cent_cpd = NaN(1,num_ev_type);
derived_spf_fwhh_cpd = NaN(1,num_ev_type);
derived_spf_fwhh_oct = NaN(1,num_ev_type);
derived_ori_fwhh_deg = NaN(1,num_ev_type);
for tp=2:num_ev_type
    derived_spf_cent_cpd(tp) = opt_values{tp}(3)/width_imgPadded_deg; % cpd
    
    temp_sigmaX_cpd = opt_values{tp}(2)/width_imgPadded_deg; % cyc/img -> cpd
    temp_hwhhX  = NormalDistrib_StdDev_to_FWHH(temp_sigmaX_cpd)/2;
    derived_spf_fwhh_cpd(tp) = temp_hwhhX*2; % cpd
    derived_spf_fwhh_oct(tp) = log2(derived_spf_cent_cpd(tp)+temp_hwhhX)...
                              -log2(derived_spf_cent_cpd(tp)-temp_hwhhX); % oct
              
    temp_sigmaY = opt_values{tp}(1);
    temp_sigmaX = opt_values{tp}(2);
    
    temp_findHWHH = @(posO) ( (0.5 -  gauss1D(opt_values{tp}(3)*sin(posO/180*pi),temp_sigmaY,0)...
                                    .*gauss1D(opt_values{tp}(3)*cos(posO/180*pi),temp_sigmaX,+opt_values{tp}(3))...
                                   -  gauss1D(opt_values{tp}(3)*sin(posO/180*pi),temp_sigmaY,0)...
                                    .*gauss1D(opt_values{tp}(3)*cos(posO/180*pi),temp_sigmaX,-opt_values{tp}(3)))^2 );
                           
    temp_value = 0; opt_value = fminsearch(temp_findHWHH,temp_value);
    
    derived_ori_fwhh_deg(tp) = 2*abs(opt_value);
end


%% Set 'plot' function
retval_mainResults.name = 'Reverse correlation method';
retval_mainResults.data_RField = receptiveFields; % [y ,x ,neTyp]
retval_mainResults.data_Fourier = fourierRF(regionY,regionX,:); % [y',x',neTyp]
retval_mainResults.data_fitted  = fittedFourierRF(regionY,regionX,:);% [y',x',neTyp]
retval_mainResults.bandwidth_spf_cpd = derived_spf_fwhh_cpd;
retval_mainResults.bandwidth_spf_oct = derived_spf_fwhh_oct;
retval_mainResults.bandwidth_ori_deg = derived_ori_fwhh_deg;
retval_mainResults.cent_spf_cpd = derived_spf_cent_cpd;
retval_mainResults.coordinatesStyle = stS.coordinatesStyle; % 'ij' or 'xy'
retval_mainResults.num_ev_type = num_ev_type;

cm1=[0:255; 0:255; 255*ones(1,256)];
cm2=[255*ones(1,256); 255:-1:0; 255:-1:0];
cm3=[cm1';cm2'];
cm3=cm3/255;
retval_mainResults.colormap = cm3;

%- Set 'plot' function
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'for tp=2:retval_mainResults.num_ev_type;'...
                                  'subplot(3,(retval_mainResults.num_ev_type-1),tp-1,''align'');'...
                                  'imagesc(retval_mainResults.data_RField(:,:,tp));'...
                                  'set(gca,''xtick'',[],''ytick'',[]); axis(retval_mainResults.coordinatesStyle); axis(''image''); colorbar;'...
                                  'subplot(3,(retval_mainResults.num_ev_type-1),tp-1+(retval_mainResults.num_ev_type-1),''align'');'...
                                  'imagesc(retval_mainResults.data_Fourier(:,:,tp));'...
                                  'set(gca,''xtick'',[],''ytick'',[]); axis(retval_mainResults.coordinatesStyle); axis(''image''); colorbar;'...
                                  'subplot(3,(retval_mainResults.num_ev_type-1),tp-1+(retval_mainResults.num_ev_type-1)*2,''align'');'...
                                  'imagesc(retval_mainResults.data_fitted(:,:,tp));'...
                                  'set(gca,''xtick'',[],''ytick'',[]); axis(retval_mainResults.coordinatesStyle); axis(''image''); colorbar;'...
                                'end;'...
                                ...
                                'fprintf(''Frequency bandwidth (cpd): '');   disp(retval_mainResults.bandwidth_spf_cpd);'...
                                'fprintf(''Frequency bandwidth (oct): '');   disp(retval_mainResults.bandwidth_spf_oct);'...
                                'fprintf(''Frequency center    (cpd): '');   disp(retval_mainResults.cent_spf_cpd);'...
                                'fprintf(''Orientation bandwidth (deg): ''); disp(retval_mainResults.bandwidth_ori_deg);'...
                            ];
retval_mainResults.plot = @(tp) eval( retval_mainResults.command );

%- Set 'comment'
retval_mainResults.comment = sprintf( 'Plot the results of ''Explore_Simple_ReverseCorrelation'' by calling ''retval_mainResults.plot()''.' );

end %%% of file
