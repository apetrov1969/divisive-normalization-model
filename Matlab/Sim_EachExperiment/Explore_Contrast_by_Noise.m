function [retval_mainResults,...
          retval_allResults] = Explore_Contrast_by_Noise(M, specs_simulation)

%Explore_Contrast_by_Noise - Contrast sensitivity (Noise)
%
% Visual stimuli:
%  Luminance gratings with various contrasts + white_noise with various contrasts
%
% References:
%  Carandini & HeegerMovshon (1997)
%  Squatrito, Trotter, & Poggio (1990)
%
% References for "Stochastic Resonance":
%  Moss, Ward, & Sannita (2004)
%  Sasaki, et al. (2006, 2008)
%  Simonotto et al. (1997)
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2);
maxGridR_deg = max(abs(gridX(:)));   %#ok

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
specs_images.radius_deg = specs_simulation.diameter_deg/2;

specs_images.orient_deg = neuron_orient_deg;
specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.wavelength_deg = 1/specs_images.spFreq_cpd; % 1/2=0.5
specs_images.size_noise_block_deg = specs_images.wavelength_deg/4; % p.8622, CarandiniHeegerMovshon_1997_JNeurosci

specs_images.condition_contrast2 = [0,6,12,25,50]./100;
specs_images.condition_contrast1_log10 = 0:0.1:log10(100); % contrast can be 50% at maximum because contrast2 can be also 50%
specs_images.condition_contrast1 = 10.^specs_images.condition_contrast1_log10 ./100; % for smooth plotting curves
specs_images.condition_contrast1 = [specs_images.condition_contrast1, 1-specs_images.condition_contrast2];
specs_images.condition_contrast1 = sort(unique(specs_images.condition_contrast1));

numConditions = [size(specs_images.condition_contrast1,2),...
                 size(specs_images.condition_contrast2,2),2,10]; % [imCon1,imCon2,imgBlock,imRand]

%- Generate random noise block frames
size_noise_block_pix = specs_images.size_noise_block_deg/stS.degPerPixel;
size_noise_frame_block = floor(stS.imageSize_pix./size_noise_block_pix);
frames_noise = NaN([stS.imageSize_pix, numConditions(4)]);
for cf=1:numConditions(4) % framses of random noise
    tempImg = rand(size_noise_frame_block);
    tempImg(tempImg<0.5)=-stS.rangeLum/2;
    tempImg(tempImg>0.0)=+stS.rangeLum/2;
    frames_noise(:,:,cf) = imresize(tempImg,[stS.imageSize_pix],'nearest');
end

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imCon1,imCon2,imgBlock,imRand]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.orientation_deg = specs_images.orient_deg;
specs_grating.frequency = specs_images.spFreq_cpd;
for c1=1:numConditions(1) % contrast1 (abscissa)
    
    % Session 1: Signal * Noise
    specs_grating.amplitude = stS.rangeLum/2 .* specs_images.condition_contrast1(c1);
    baseGrating = Grating2D(specs_grating,gridX,gridY);
    for c2=1:numConditions(2) % contrast2
        for cf=1:numConditions(4) % framses of random noise
            tempImg = baseGrating + stS.bgLum + squeeze(frames_noise(:,:,cf))*specs_images.condition_contrast2(c2);
            tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
            setStimuli(:,:,c1,c2,1,cf)=tempImg;
            if specs_images.condition_contrast1(c1) + specs_images.condition_contrast2(c2) > 1.001
                setStimuli(:,:,c1,c2,1,cf)=0;
            end
        end
    end
    
    % Session 2: Noise * Signal
    for c2=1:numConditions(2) % contrast2
        specs_grating.amplitude = stS.rangeLum/2 .* specs_images.condition_contrast2(c2);
        baseGrating = Grating2D(specs_grating,gridX,gridY);
        for cf=1:numConditions(4) % framses of random noise
            tempImg = baseGrating + stS.bgLum + squeeze(frames_noise(:,:,cf))*specs_images.condition_contrast1(c1);
            tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
            setStimuli(:,:,c1,c2,2,cf)=tempImg;
            if specs_images.condition_contrast1(c1) + specs_images.condition_contrast2(c2) > 1.001
                setStimuli(:,:,c1,c2,2,cf)=0;
            end
        end
    end
end
setStimuli = reshape(setStimuli,[stS.imageSize_pix, numConditions(1)*numConditions(2)*numConditions(3)*numConditions(4)]); % [y,x,img]


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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2),numConditions(3),numConditions(4)]); % [neTyp,imCon1,imCon2,imgBlock,imRand]
results = mean(results,5); % [neTyp,imCon1,imCon2,imgBlock,1]
results_signal = squeeze(results(:,:,:,1,:)); % [neTyp,imCon1,imCon2]
results_noise  = squeeze(results(:,:,:,2,:)); % [neTyp,imCon1,imCon2]

%- Remove contrast > 100%
for c1=1:size(results_signal,2) % contrast1
    for c2=1:size(results_signal,3) % contrast2
        if(specs_images.condition_contrast1(c1) + ...
           specs_images.condition_contrast2(c2) > 1.0001)
            results_signal(:,c1,c2) = NaN;
            results_noise (:,c1,c2) = NaN;
        end
    end % for c2 % contrast2
end % for c1 % contrast1

%- Set 'retval_mainResults'
retval_mainResults.name = 'Contrast x Noise';
retval_mainResults.data1 = results_signal;
retval_mainResults.data2 = results_noise;
retval_mainResults.xval = specs_images.condition_contrast1;
retval_mainResults.xlabel1 = 'Grating Contrast';
retval_mainResults.xlabel2 = 'Noise Contrast';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items1 = cellstr(num2str(specs_images.condition_contrast2','%.2f'));
retval_mainResults.legend_items2 = cellstr(num2str(specs_images.condition_contrast2','%.2f'));
retval_mainResults.legend_loc   = 'SouthEast';

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'subplot(1, 2, 1, ''align'');'...
                                'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.data1(temp_type,:,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel1);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items1, ''Location'',retval_mainResults.legend_loc);'...
                                ...
                                'subplot(1, 2, 2, ''align'');'...
                                'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.data2(temp_type,:,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel2);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items2, ''Location'',retval_mainResults.legend_loc)'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Contrast_by_Noise'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );
                                   
end %%% of file
