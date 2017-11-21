function [retval_mainResults,...
          retval_allResults] = Explore_XOrient2(M, specs_simulation)

      
      
      
      
      %Base and mask spfs can be controlled.
      
      
      
      
      
      
      
%Explore_XOrient - X orientation tuning
%
% Visual stimuli:
%  Two gratings superimposed with various orientations
%
% References:
%  De Angelis, Robson, Ohzawa, & Freeman (1992)
%  Morrone et al. (1982)
%  Bonds_1989
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2);
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
specs_images.radius_deg  = specs_simulation.diameter_deg/2;
specs_images.base_contrast = specs_simulation.contrast1_pcent/100;
specs_images.mask_contrast  = specs_simulation.contrast2_pcent/100;

specs_images.base_spFreq_cpd = specs_simulation.base_spFreq_cpd;
specs_images.mask_spFreq_cpd = specs_simulation.mask_spFreq_cpd;

specs_images.base_orient_deg = neuron_orient_deg;

specs_images.condition_var_orient_deg = -90:5:+90;
numConditions=[size(specs_images.condition_var_orient_deg,2),2]; % [imOri,single/x]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imOri+1,single/x]
specs_baseGrating.type = 'Cos';
specs_baseGrating.phase_deg = 0;
specs_baseGrating.frequency = specs_images.base_spFreq_cpd;
specs_baseGrating.amplitude = stS.rangeLum/2 .* specs_images.base_contrast;
specs_baseGrating.orientation_deg = specs_images.base_orient_deg;
baseGrating = Grating2D(specs_baseGrating,gridX,gridY); 

specs_varGrating.type = 'Cos';
specs_varGrating.phase_deg = 90;
for c=1:numConditions(1) % orientation
    specs_varGrating.orientation_deg = specs_images.condition_var_orient_deg(c);
    
    %-- Single grating
    specs_varGrating.amplitude = stS.rangeLum/2 .* specs_images.base_contrast;
    specs_varGrating.frequency = specs_images.base_spFreq_cpd;
    varGrating = Grating2D(specs_varGrating,gridX,gridY);
	tempImg = varGrating + stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
    setStimuli(:,:,c,1) = tempImg;
    
    %-- X grating
    specs_varGrating.amplitude = stS.rangeLum/2 .* specs_images.mask_contrast;
    specs_varGrating.frequency = specs_images.mask_spFreq_cpd;
    varGrating = Grating2D(specs_varGrating,gridX,gridY);
	tempImg = baseGrating + varGrating + stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
    setStimuli(:,:,c,2) = tempImg;
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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imOri,single/x]

%- Set 'retval_mainResults'
retval_mainResults.name = 'X-orientation suppression';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_var_orient_deg;
retval_mainResults.xlabel = 'Orientation (deg)';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'Single';'Cross'};
retval_mainResults.legend_loc   = 'NorthWest';
retval_mainResults.xtick = -90:45:90;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'plot(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:,:)),''-'');'...
                                'set(gca,''xtick'',retval_mainResults.xtick);'...
                                'title(retval_mainResults.name);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_XOrient'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
