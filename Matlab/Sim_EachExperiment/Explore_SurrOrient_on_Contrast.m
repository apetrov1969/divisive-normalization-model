function [retval_mainResults,...
          retval_allResults] = Explore_SurrOrient_on_Contrast(M, specs_simulation)

%Explore_SurrOrient_on_Contrast - Contrast sensitivity (surround orientation)
%
% Visual stimuli:
%  Luminance gratings with various contrasts in center and
%                     with various orientations in surround
%
% References:
%  Cavanagh, Bair, & Movshon (2002)
%  Levitt & Lund (1997): Facilitation
%  Polat et al.  (1998): Facilitation
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2);
maxGridR_deg = max(abs(gridX(:)));

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
specs_images.contrast2 = specs_simulation.surr_contrast_pcent/100;
specs_images.radius_Cent_deg = specs_simulation.center_diameter_deg/2;
specs_images.radius_inner_Surr_deg = specs_simulation.surround_inner_diameter_out_deg/2;
specs_images.radius_outer_Surr_deg = specs_simulation.surround_outer_diameter_out_deg/2;

specs_images.orient_deg = neuron_orient_deg;
specs_images.spFreq_cpd = neuron_spFreq_cpd;

specs_images.condition_contrast1_log10 = 0:0.1:2;
specs_images.condition_contrast1 = 10.^specs_images.condition_contrast1_log10 ./100;
numConditions = [size(specs_images.condition_contrast1,2),3]; % [imCon,imTyp]
[dummy,idx_cont100_img] = min(abs(specs_images.condition_contrast1-1.0));  %#ok
[dummy,idx_cont050_img] = min(abs(specs_images.condition_contrast1-0.5));  %#ok
[dummy,idx_cont012_img] = min(abs(specs_images.condition_contrast1-0.12)); %#ok

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imOri,imType]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.frequency = specs_images.spFreq_cpd;

% Orthogonal surround grating (Horizontal)
specs_grating.orientation_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.contrast2;
orthGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;

% Parallel surround grating (Vertical)
specs_grating.orientation_deg = 90;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.contrast2;
paraGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;

specs_grating.orientation_deg = specs_images.orient_deg;
for c1=1:numConditions(1) % contrast1: center
    specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast1(c1);
    centGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
        
    % Type1: center grating only
    tempImg=centGrating;
    tempImg(gridR_deg>(specs_images.radius_Cent_deg))=stS.bgLum;
    setStimuli(:,:,c1,1) = tempImg;
    
    % Type2: orthogonal surround grating    
    tempImg=orthGrating;
    tempImg(gridR_deg>specs_images.radius_outer_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_inner_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_Cent_deg)=centGrating(gridR_deg<(specs_images.radius_Cent_deg));
    setStimuli(:,:,c1,2) = tempImg;
    
    % Type3: parallel surround grating
    tempImg=paraGrating;
    tempImg(gridR_deg>specs_images.radius_outer_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_inner_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_Cent_deg)=centGrating(gridR_deg<(specs_images.radius_Cent_deg));
    setStimuli(:,:,c1,3) = tempImg;
end
setStimuli = reshape(setStimuli,[stS.imageSize_pix, numConditions(1)*numConditions(2)]); % [y,x,img]


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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imCon,imTyp]

%- Compute Suppression Ratio
suppressionRatio=NaN(num_ev_type,3,3); % Contrast (100/50/12) x Surround grating (Orth/Para)
for tp=1:num_ev_type
    suppressionRatio(tp,1,1)=100; % Contrast 100%
    suppressionRatio(tp,1,2)=results(tp,idx_cont100_img,2)/results(tp,idx_cont100_img,1);
    suppressionRatio(tp,1,3)=results(tp,idx_cont100_img,3)/results(tp,idx_cont100_img,1);

    suppressionRatio(tp,2,1)=50; % Contrast 50%
    suppressionRatio(tp,2,2)=results(tp,idx_cont050_img,2)/results(tp,idx_cont050_img,1);
    suppressionRatio(tp,2,3)=results(tp,idx_cont050_img,3)/results(tp,idx_cont050_img,1);

    suppressionRatio(tp,3,1)=12; % Contrast 12%
    suppressionRatio(tp,3,2)=results(tp,idx_cont012_img,2)/results(tp,idx_cont012_img,1);
    suppressionRatio(tp,3,3)=results(tp,idx_cont012_img,3)/results(tp,idx_cont012_img,1);
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Surround suppression (Contrast x Orient)';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_contrast1;
retval_mainResults.xlabel = 'Center Contrast (%)';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'Cent';'Orth';'Para'};
retval_mainResults.legend_loc   = 'NorthWest';
retval_mainResults.suppressionRatio = suppressionRatio;


%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_SurrOrient_on_Contrast'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
