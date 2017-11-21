function [retval_mainResults,...
          retval_allResults] = Explore_Orient_by_Size(M, specs_simulation)

%Explore_Orient_by_Size - Orientation tuning (size)
%
% Visual stimuli:
%  Luminance gratings with various orientations and size
%
% Refrerences
%  Okamoto, Naito, Sadakane, Osaki, & Sato (2009)
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
specs_images.luminance_contrast = specs_simulation.contrast_pcent/100;
specs_images.base_radius_deg = specs_simulation.diameter_deg/2;

if(specs_images.base_radius_deg*4>maxGridR_deg)
    warning('M.stim_spec.imageSize_pix is too small for ''Explore_Orient_by_Size''');
end

specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.condition_radius_deg = [specs_images.base_radius_deg*1,...
                                     specs_images.base_radius_deg*2,...
                                     specs_images.base_radius_deg*4];
specs_images.condition_diameter_deg = specs_images.condition_radius_deg*2;
specs_images.orient_step_deg = 2;
specs_images.condition_orient_deg = -90:specs_images.orient_step_deg:90;
numConditions = [size(specs_images.condition_orient_deg,2),...
                 size(specs_images.condition_diameter_deg,2)]; % [imOri,imSize]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions(1), numConditions(2)]); % [x,y,imOri,imSize]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast;
specs_grating.frequency = specs_images.spFreq_cpd;
for c1=1:numConditions(1) % orient
    for c2=1:numConditions(2) % size (diameter)
        specs_grating.orientation_deg = specs_images.condition_orient_deg(c1);
        tempImg = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
        tempImg(gridR_deg>specs_images.condition_radius_deg(c2))=stS.bgLum;
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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imOri,imSize]

bandwidth  =NaN(num_ev_type,numConditions(2));
fullH_idx  =NaN(num_ev_type,1,numConditions(2));
halfH_idx  =NaN(num_ev_type,2,numConditions(2));
fullH_yval =NaN(num_ev_type,1,numConditions(2));
halfH_yval =NaN(num_ev_type,1,numConditions(2));
fullH_xval =NaN(num_ev_type,1,numConditions(2));
halfH_xval =NaN(num_ev_type,2,numConditions(2));

for c2=1:numConditions(2) % size (diameter)
    for tp=1:num_ev_type
        [fullH_yval(tp,1,c2),fullH_idx(tp,1,c2),...
         halfH_yval(tp,1,c2),halfH_idx(tp,:,c2)]=DeriveFWHH(results(tp,:,c2));

        fullH_xval(tp,1,c2) = ValueFloatIdx(specs_images.condition_orient_deg,fullH_idx(tp,1,c2));
        halfH_xval(tp,1,c2) = ValueFloatIdx(specs_images.condition_orient_deg,halfH_idx(tp,1,c2));
        halfH_xval(tp,2,c2) = ValueFloatIdx(specs_images.condition_orient_deg,halfH_idx(tp,2,c2));

        bandwidth(tp,c2) = halfH_xval(tp,2,c2)-halfH_xval(tp,1,c2);
    end
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Orientation x Size';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_orient_deg;
retval_mainResults.xlabel = 'Orientation';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = cellstr(num2str(specs_images.condition_diameter_deg','%.1f�'));
retval_mainResults.legend_loc   = 'NorthWest';
retval_mainResults.xtick = -90:45:90;
retval_mainResults.bandwidth = bandwidth;

retval_mainResults.marking_y = NaN(num_ev_type,3,numConditions(2));
retval_mainResults.marking_y(:,1,:) = fullH_yval;
retval_mainResults.marking_y(:,2,:) = halfH_yval;
retval_mainResults.marking_y(:,3,:) = halfH_yval;

retval_mainResults.marking_x = NaN(num_ev_type,3,numConditions(2));
retval_mainResults.marking_x(:,1,:)   = fullH_xval;
retval_mainResults.marking_x(:,2:3,:) = halfH_xval;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'plot(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'set(gca,''xtick'',retval_mainResults.xtick);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                                'hold on;'...
                                'plot(squeeze(retval_mainResults.marking_x(temp_type,:,:)),  squeeze(retval_mainResults.marking_y(temp_type,:,:)),  ''o'');'...
                                'hold off;'...
                                'retval_mainResults.bandwidth(temp_type,:)'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Orient_by_Size'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
