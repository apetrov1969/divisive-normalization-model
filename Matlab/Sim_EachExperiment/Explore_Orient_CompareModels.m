function [retval_mainResults,...
          retval_allResults] = Explore_Orient_CompareModels(M, specs_simulation)

%Explore_Orient - Orientation tuning of models
%
% Visual stimuli:
%  Luminance gratings with various orientations
%
% Refrerences
%  Okamoto, Naito, Sadakane, Osaki, & Sato (2009)
%  Rose & Blakemore (1974)
%  Watkins & Berkley (1974)
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
%num_ev_type = evS.num_cellType; % [1 Complex cell, 4 Simple cells]


%% Set up simulation
%- Get indices of neurons for analysing the results
neuron_orient_deg = specs_simulation.neuron_orient_deg;
neuron_spFreq_cpd = specs_simulation.neuron_spFreq_cpd;
idx_neuron_orient = FindClosestIdx(evS.domain_orient_deg,neuron_orient_deg);
idx_neuron_spFreq = FindClosestIdx(evS.domain_spFreq_cpd,neuron_spFreq_cpd);
idx_neuron_rfLoc = 1;

%- Set specs of stimuli
specs_images.luminance_contrast = specs_simulation.contrast_pcent/100;
specs_images.radius_deg = specs_simulation.diameter_deg/2;

specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.orient_step_deg = 2;
specs_images.condition_orient_deg = -90:specs_images.orient_step_deg:90;
numConditions = size(specs_images.condition_orient_deg,2);

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imOri]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast;
specs_grating.frequency = specs_images.spFreq_cpd;
for c1=1:numConditions % orient
    specs_grating.orientation_deg = specs_images.condition_orient_deg(c1);
    tempImg = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
    tempImg(gridR_deg > specs_images.radius_deg) = stS.bgLum;
    setStimuli(:,:,c1) = tempImg;
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
idx_neuron_type = 1; % complex
results_dnm = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,idx_neuron_type,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [1,imOri]
results_Stim = squeeze(resp_StimDrv(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,idx_neuron_type,:)); % [1,imOri]

results_dnm  = reshape(results_dnm ,[1,numConditions]);
results_Stim = reshape(results_Stim,[1,numConditions]);

exponent     = M.EarlyVis_spec.divNorm_exponentNum; % 2.0
%semisatConst = M.EarlyVis_spec.divNorm_semisaturConst; % 0.1000
height = max(results_dnm(:));
results_Power = results_Stim.^exponent;
results = [height.*results_Stim; height.*results_Power; results_dnm];
num_models = 3;


bandwidth  =NaN(num_models,1);
fullH_idx  =NaN(num_models,1);
halfH_idx  =NaN(num_models,2);
fullH_yval =NaN(num_models,1);
halfH_yval =NaN(num_models,1);
fullH_xval =NaN(num_models,1);
halfH_xval =NaN(num_models,2);

for tp=1:num_models
    [fullH_yval(tp,1),fullH_idx(tp,1),...
     halfH_yval(tp,1),halfH_idx(tp,:)]=DeriveFWHH(results(tp,:));
 
    fullH_xval(tp,1) = ValueFloatIdx(specs_images.condition_orient_deg,fullH_idx(tp,1));
    halfH_xval(tp,1) = ValueFloatIdx(specs_images.condition_orient_deg,halfH_idx(tp,1));
    halfH_xval(tp,2) = ValueFloatIdx(specs_images.condition_orient_deg,halfH_idx(tp,2));
    
    bandwidth(tp) = halfH_xval(tp,2)-halfH_xval(tp,1);
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Orientation';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_orient_deg;
retval_mainResults.xlabel = 'Orientation';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.xtick = -90:45:90;
retval_mainResults.legend_items = {'Linear-Rect';'Exponential';'Divisive-Normalization'};

retval_mainResults.bandwidth_deg = bandwidth;
retval_mainResults.fullH_idx = fullH_idx;
retval_mainResults.halfH_idx = halfH_idx;
retval_mainResults.fullH_xval = fullH_xval;
retval_mainResults.halfH_xval = halfH_xval;
retval_mainResults.fullH_yval = fullH_yval;
retval_mainResults.halfH_yval = halfH_yval;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'plot(retval_mainResults.xval,squeeze(retval_mainResults.data(:,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'set(gca,''xtick'',retval_mainResults.xtick);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                ...%'axis([-90,90,0,inf]);'...
                                'legend(retval_mainResults.legend_items);'...
                                ...%'hold on;'...
                                ...%'plot(retval_mainResults.fullH_xval(temp_type),  retval_mainResults.fullH_yval(temp_type),  ''o'');'...
                                ...%'plot(retval_mainResults.halfH_xval(temp_type,1),retval_mainResults.halfH_yval(temp_type,1),''o'');'...
                                ...%'plot(retval_mainResults.halfH_xval(temp_type,2),retval_mainResults.halfH_yval(temp_type,1),''o'');'...
                                ...%'hold off;'...
                                'disp(retval_mainResults.legend_items(:));'...
                                'disp(retval_mainResults.bandwidth_deg(:));'...
                            ];
retval_mainResults.plot = @() eval( retval_mainResults.command );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Orient_CompareModels'' by calling ''retval_mainResults.plot()''.\n'] );

end %%% of file
