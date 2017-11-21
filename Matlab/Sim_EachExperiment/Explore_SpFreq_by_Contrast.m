function [retval_mainResults,...
          retval_allResults] = Explore_SpFreq_by_Contrast(M, specs_simulation)

%Explore_SpFreq_by_Contrast - Spatial-frequency tuning (contrast)
%
% Visual stimuli:
%  Luminance gratings with various frequencies and contrasts
%
% References:
%  Sceniak, Hawken, & Shapley (2001)
%  Skottun, et al. (1987)
%  Albrecht & Hamilton (1982)
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
specs_images.radius_deg = specs_simulation.diameter_deg/2;

specs_images.orient_deg = neuron_orient_deg;
specs_images.condition_contrast = [0.05,0.10,0.20,0.40,0.80];
specs_images.condition_spFreq_oct = -1:0.01:+3;
specs_images.condition_spFreq_cpd = 2.^(specs_images.condition_spFreq_oct);
numConditions = [size(specs_images.condition_spFreq_cpd,2),...
                 size(specs_images.condition_contrast,2)];

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.orientation_deg = specs_images.orient_deg;
for c1=1:numConditions(1) % spFreq
    for c2=1:numConditions(2) % contrast
        specs_grating.frequency = specs_images.condition_spFreq_cpd(c1);
        specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast(c2);
        tempImg = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
        tempImg(gridR_deg > specs_images.radius_deg) = stS.bgLum;
        setStimuli(:,:,c1,c2) = tempImg;
    end
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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imSpf,imCon]

bandwidth_cpd=NaN(num_ev_type,numConditions(2));
bandwidth_oct=NaN(num_ev_type,numConditions(2));

fullH_idx  =NaN(num_ev_type,1,numConditions(2));
halfH_idx  =NaN(num_ev_type,2,numConditions(2));
fullH_yval =NaN(num_ev_type,1,numConditions(2));
halfH_yval =NaN(num_ev_type,1,numConditions(2));
fullH_xval_cpd =NaN(num_ev_type,1,numConditions(2));
halfH_xval_cpd =NaN(num_ev_type,2,numConditions(2));
fullH_xval_oct =NaN(num_ev_type,1,numConditions(2));
halfH_xval_oct =NaN(num_ev_type,2,numConditions(2));

for c2=1:numConditions(2) % contrast
    for tp=1:num_ev_type
        [fullH_yval(tp,1,c2),fullH_idx(tp,1,c2),...
         halfH_yval(tp,1,c2),halfH_idx(tp,:,c2)]=DeriveFWHH(squeeze(results(tp,:,c2)));
        
        fullH_xval_cpd(tp,1,c2) = ValueFloatIdx(specs_images.condition_spFreq_cpd,fullH_idx(tp,1,c2));
    	halfH_xval_cpd(tp,1,c2) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,1,c2));
        halfH_xval_cpd(tp,2,c2) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,2,c2));
    
        fullH_xval_oct(tp,1,c2) = ValueFloatIdx(specs_images.condition_spFreq_oct,fullH_idx(tp,1,c2));
    	halfH_xval_oct(tp,1,c2) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,1,c2));
        halfH_xval_oct(tp,2,c2) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,2,c2));
        
        bandwidth_cpd(tp,c2) = halfH_xval_cpd(tp,2,c2)-halfH_xval_cpd(tp,1,c2);
        bandwidth_oct(tp,c2) = halfH_xval_oct(tp,2,c2)-halfH_xval_oct(tp,1,c2);
    end
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Spatial frequency x Contrast';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_spFreq_oct;
retval_mainResults.xlabel = 'Spatial frequency (deg)';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = cellstr(num2str(specs_images.condition_contrast'*100,'%d%%°'));
retval_mainResults.legend_loc   = 'NorthWest';

retval_mainResults.bandwidth_cpd = bandwidth_cpd;
retval_mainResults.bandwidth_oct = bandwidth_oct;

retval_mainResults.marking_y = NaN(num_ev_type,3,numConditions(2));
retval_mainResults.marking_y(:,1,:) = fullH_yval;
retval_mainResults.marking_y(:,2,:) = halfH_yval;
retval_mainResults.marking_y(:,3,:) = halfH_yval;

retval_mainResults.marking_x = NaN(num_ev_type,3,numConditions(2));
retval_mainResults.marking_x(:,1,:)   = fullH_xval_oct;
retval_mainResults.marking_x(:,2:3,:) = halfH_xval_oct;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'plot(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                                'hold on;'...
                                'plot(squeeze(retval_mainResults.marking_x(temp_type,:,:)),squeeze(retval_mainResults.marking_y(temp_type,:,:)),''o'');'...
                                'hold off;'...
                                'disp(retval_mainResults.bandwidth_oct(temp_type,:));'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_SpFreq_by_Contrast'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
