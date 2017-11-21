function [retval_resp_centComp,...
          retval_resp_centSimp,...
          retval_resp_suppChan,...
          retval_resp_energyMap] = DMPL_EarlyVis_PoolRectifiedLinChannels(M,stimuli)


assert(M.prepare_level >= 2);
assert(M.EarlyVis_spec.prepare_level >= 3);

assert(size(stimuli,1)==M.stim_spec.imageSize_pix(1));
assert(size(stimuli,2)==M.stim_spec.imageSize_pix(2));


%% Retrieve parameters
num_chan_SpFreq = M.EarlyVis_spec.num_channel_spFreq; % chanSpf
num_orient = M.EarlyVis_spec.num_orient; % neurOri/chanOri
num_spFreq = M.EarlyVis_spec.num_spFreq; % neurSpf
num_rfLoc  = M.EarlyVis_spec.num_rfLoc;  % neurRfl
num_images = size(stimuli,3);

size_imgPadded = M.EarlyVis_filters.size_imgPadded;
size_padRegion = M.EarlyVis_filters.size_leftBtmPad;

imageSize_pix = M.stim_spec.imageSize_pix;
num_pixels = prod(imageSize_pix);

divNorm_exponentDen = M.EarlyVis_spec.divNorm_exponentDen;


%% Prepare to pad stimuli with zeros to avoid FFT-related convolution contamination
regionY = (1:imageSize_pix(1))+size_padRegion(1); %y: 1st index
regionX = (1:imageSize_pix(2))+size_padRegion(2); %x: 2nd index
paddedStimulus = zeros(size_imgPadded);


%% Rehape matrices for image processing
filters_pool_surr_yx = M.EarlyVis_poolKernels.filters_pool_surr_yx;
filters_pool_surr_yx = reshape(filters_pool_surr_yx,[num_rfLoc*num_spFreq,num_pixels]); % [neurRfl*neurSpf,yx]
clrf_cos_filter = reshape(M.EarlyVis_filters.filters_clrf_CosImg,[num_pixels,num_orient*num_spFreq*num_rfLoc]); % [yx,neurOri*neurSpf*neurRfl]
clrf_sin_filter = reshape(M.EarlyVis_filters.filters_clrf_SinImg,[num_pixels,num_orient*num_spFreq*num_rfLoc]); % [yx,neurOri*neurSpf*neurRfl]


%% Apply the image filters (using Convolution theorem), rectify, and pool

%- Convert luminance to contrast with mean luminance = bgLum:  [minLum maxLum] -> [-0.5 +0.5]
stimContrast = (stimuli-M.stim_spec.bgLum) ./ M.stim_spec.rangeLum;

%- Allocate memory
retval_resp_suppChan  = NaN(num_orient,num_chan_SpFreq,num_spFreq,num_rfLoc,num_images); % [chanOri, chanSpf, neurSpf, neurRfl, img]
retval_resp_centComp  = NaN(num_orient,num_spFreq,num_rfLoc,num_images);                 % [neurOri, neurSpf, neurRfl, img]
retval_resp_centSimp  = NaN(num_orient,num_spFreq,num_rfLoc,4,num_images);               % [neurOri, neurSpf, neurRfl, neurPhs, img]
if(nargout>=4)
    retval_resp_energyMap = NaN(imageSize_pix(1),imageSize_pix(2),num_orient,num_chan_SpFreq,num_images); % [y,x,chanOri,chanSpf, img]
end

for s = 1:num_images
    %% Surround suppressive signal: Convolution
    % Use Convolution theorem:
    %   Convolution(stim,filter)      == ifft(fft(stim) .* fft(filter))
    %   CrossCorrelation(stim,filter) == ifft(fft(stim) .* conj(fft(filter)))
    paddedStimulus(:) = 0 ; % fill with fresh zeros just in case
    paddedStimulus(regionY,regionX) = stimContrast(:,:,s);
    paddedStimulusFFT = fft2(paddedStimulus);
    
    energyMap4D = zeros(imageSize_pix(1),imageSize_pix(2),num_orient,num_chan_SpFreq);
    if(strcmp(M.EarlyVis_spec.phaseInvariance_method,'quadrature'))
        for o = 1:num_orient
            for f = 1:num_chan_SpFreq
                stimulusFilteredFFTCos = paddedStimulusFFT .* M.EarlyVis_filters.filters_Surr_CosFFT{o,f};
                mapCos = real(ifft2(stimulusFilteredFFTCos));
                stimulusFilteredFFTSin = paddedStimulusFFT .* M.EarlyVis_filters.filters_Surr_SinFFT{o,f};
                mapSin = real(ifft2(stimulusFilteredFFTSin));
                mapAbs = sqrt(mapCos.^2 + mapSin.^2);
                energyMap4D(:,:,o,f) = mapAbs(regionY,regionX); % [y,x,chanOri,chanSpf]
            end % for f = 1:num_chan_SpFreq
        end % for o = 1:num_orient
    elseif(strcmp(M.EarlyVis_spec.phaseInvariance_method,'cos-only'))
        for o = 1:num_orient
            for f = 1:num_chan_SpFreq
                stimulusFilteredFFTCos = paddedStimulusFFT .* M.EarlyVis_filters.filters_Surr_CosFFT{o,f};
                mapCos = real(ifft2(stimulusFilteredFFTCos));
                mapAbs = abs(mapCos);
                energyMap4D(:,:,o,f) = mapAbs(regionY,regionX); % [y,x,chanOri,chanSpf]
            end % for f = 1:num_chan_SpFreq
        end % for o = 1:num_orient
    elseif(strcmp(M.EarlyVis_spec.phaseInvariance_method,'sin-only'))
        for o = 1:num_orient
            for f = 1:num_chan_SpFreq
                stimulusFilteredFFTSin = paddedStimulusFFT .* M.EarlyVis_filters.filters_Surr_SinFFT{o,f};
                mapSin = real(ifft2(stimulusFilteredFFTSin));
                mapAbs = abs(mapSin);
                energyMap4D(:,:,o,f) = mapAbs(regionY,regionX); % [y,x,chanOri,chanSpf]
            end % for f = 1:num_chan_SpFreq
        end % for o = 1:num_orient
    else
        error('Invalid M.EarlyVis_spec.phaseInvariance_method.');
    end

    if(nargout>=4)
        retval_resp_energyMap(:,:,:,:,s) = energyMap4D; % [y,x,chanOri,chanSpf, img]
    end
    energyMap4D = energyMap4D.^divNorm_exponentDen; % [y,x,chanOri,chanSpf]
	energyMap4D = reshape(energyMap4D,[num_pixels,num_orient*num_chan_SpFreq]);                     % [yx,chanOri*chanSpf]
    pooledEnergySurr = filters_pool_surr_yx * energyMap4D;                                           % [neurRfl*neurSpf,chanOri*chanSpf]
    pooledEnergySurr = reshape(pooledEnergySurr, [num_rfLoc,num_spFreq,num_orient,num_chan_SpFreq]); % [neurRfl,neurSpf,chanOri,chanSpf]
    pooledEnergySurr = permute(pooledEnergySurr, [3,4,2,1]);                                         % [chanOri,chanSpf,neurSpf,neurRfl]
    retval_resp_suppChan(:,:,:,:,s) = pooledEnergySurr;                                              % [chanOri,chanSpf,neurSpf,neurRfl,img]
    
    
    %% Simple cell: Dot-product
    stimVector = stimContrast(:,:,s);
    stimVector = stimVector(:)'; % [1,yx]
    
	simple_phaseCos = stimVector*clrf_cos_filter; % [1,neurOri*neurSpf*neurRfl]
    simple_phaseSin = stimVector*clrf_sin_filter; % [1,neurOri*neurSpf*neurRfl]
    simple_phaseCos = reshape(simple_phaseCos,[num_orient,num_spFreq,num_rfLoc]); % [neurOri,neurSpf,neurRfl]
    simple_phaseSin = reshape(simple_phaseSin,[num_orient,num_spFreq,num_rfLoc]); % [neurOri,neurSpf,neurRfl]
    
    retval_resp_centSimp(:,:,:,1,s) = +simple_phaseCos; % phase=000; % [neurOri,neurSpf,neurRfl,neurPhs,img]
    retval_resp_centSimp(:,:,:,2,s) = +simple_phaseSin; % phase=090; % [neurOri,neurSpf,neurRfl,neurPhs,img]
    retval_resp_centSimp(:,:,:,3,s) = -simple_phaseCos; % phase=180; % [neurOri,neurSpf,neurRfl,neurPhs,img]
    retval_resp_centSimp(:,:,:,4,s) = -simple_phaseSin; % phase=270; % [neurOri,neurSpf,neurRfl,neurPhs,img]
    
    
    %% Complex cell: Based on Simple cell
    retval_resp_centComp(:,:,:,s) = sqrt(simple_phaseCos.^2 + simple_phaseSin.^2); % [neurOri,neurSpf,neurRfl,img]
    
end % s = 1:num_images

%%% Return retval_PoolRectifiedLinChannels
end %%% of file