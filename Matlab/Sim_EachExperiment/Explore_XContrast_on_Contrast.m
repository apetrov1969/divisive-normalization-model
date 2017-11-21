function [retval_mainResults,...
          retval_allResults] = Explore_XContrast_on_Contrast(M, specs_simulation)

%Explore_XContrast_on_Contrast - Contrast sensitivity (surround contrast)
%
% Visual stimuli:
%  Two gratings superimposed with various contrasts
%
% References:
%  Carandini, Heeger & Movshon (1997)
%  Carandini (2004)
%  Freeman, Durand, Kiper, & Carandini (2002)
%  Morrone et al. (1982)
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2);
maxGridR_deg = max(abs(gridX(:))); %#ok

num_ev_ori  = evS.num_orient;
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
specs_images.radius_deg  = specs_simulation.diameter_deg/2;
specs_images.orient1_deg = neuron_orient_deg;

specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.orient2_deg = specs_images.orient1_deg+90;
[dummy,idx_neuron_ortho] = FindClosest(evS.domain_orient_deg,specs_images.orient2_deg); %#ok

specs_images.condition_contrast1 = 10.^(0:0.1:log10(50))./100;
specs_images.condition_contrast2 = [0,6,12,25,50]./100;
specs_images.condition_contrast1 = [specs_images.condition_contrast1, 1-specs_images.condition_contrast2];
specs_images.condition_contrast1 = sort(unique(specs_images.condition_contrast1));

numConditions=[size(specs_images.condition_contrast1,2),...
               size(specs_images.condition_contrast2,2)]; % [imCon1,imCon2]

% Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imCon1,imCon2]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.frequency = specs_images.spFreq_cpd;

for c1=1:numConditions(1) % contrast 1
    specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast1(c1);
    specs_grating.orientation_deg = specs_images.orient1_deg;
    grating1 = Grating2D(specs_grating,gridX,gridY);
    for c2=1:numConditions(2) % contrast 2
        specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast2(c2);
        specs_grating.orientation_deg = specs_images.orient2_deg;
        grating2 = Grating2D(specs_grating,gridX,gridY);
        tempImg = grating1+grating2+stS.bgLum;
        tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
        setStimuli(:,:,c1,c2) = tempImg;
    end % for c2=1:numConditions(2) % contrast 2
end % for c1=1:numConditions(1) % contrast 1
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
results = squeeze(resp_DivNorm(:,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neOri,neTyp,img]
results = reshape(results,[num_ev_ori,num_ev_type,numConditions(1), numConditions(2)]); % [neOri,neTyp,imCon1,imCon2]

%- Remove contrast > 100%
for c1=1:size(results,3) % contrast1
    for c2=1:size(results,4) % contrast2
        if(specs_images.condition_contrast1(c1) + ...
           specs_images.condition_contrast2(c2) > 1.0001)
            results(:,:,c1,c2) = NaN;
        end
    end % for c2 % contrast2
end % for c1 % contrast1

%- Set 'retval_mainResults'
retval_mainResults.name{1} = 'X-Gratings: Effect of Base contrast';
retval_mainResults.name{2} = 'X-Gratings: Effect of Mask contrast';
retval_mainResults.data = results; % [neOri,neTyp,imCon1,imCon2]
retval_mainResults.xval = specs_images.condition_contrast1;
retval_mainResults.xlabel{1} = 'Base contrast';
retval_mainResults.xlabel{2} = 'Mask contrast';
retval_mainResults.ylabel = 'Firing Rate';

retval_mainResults.neuron_idx  = [idx_neuron_orient,idx_neuron_ortho];
retval_mainResults.legend_items = cellstr(num2str(specs_images.condition_contrast2'*100,'%d%%'));
retval_mainResults.legend_loc   = 'NorthWest';

%- Set 'plot' function
% (1) Effect of Base contrast
%      Figures 2 and 10 in FreemanDurandKiperCarandini_2002_Neuron
%      Figure 10 in Carandini, Heeger & Movshon (1997)
%
% (2) Effect of Mask contrast
%      Figure 2  in FreemanDurandKiperCarandini_2002_Neuron
%      Figure 10 & 12 in Carandini, Heeger & Movshon (1997)
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
fg=1; %#ok
retval_mainResults.command = [
                                'fg=1;'...
                                'figure(''Name'',retval_mainResults.name{fg});'...
                                'semilogx(retval_mainResults.xval, squeeze(retval_mainResults.data(retval_mainResults.neuron_idx(fg),temp_type,:,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel{fg});'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                                'fg=2;'...
                                'figure(''Name'',retval_mainResults.name{fg});'...
                                'semilogx(retval_mainResults.xval, squeeze(retval_mainResults.data(retval_mainResults.neuron_idx(fg),temp_type,:,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel{fg});'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_XContrast_on_Contrast'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file