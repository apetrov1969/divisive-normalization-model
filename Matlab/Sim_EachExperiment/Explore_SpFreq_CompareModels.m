function [retval_mainResults,...
          retval_allResults] = Explore_SpFreq_CompareModels(M, specs_simulation)

%Explore_SpFreq - Spatial-frequency tuning of models
%
% Visual stimuli:
%  Luminance gratings with various spatial-frequencies
%
% Refrerences
%  Movshon, Thompson, & Tolhurst (1978)
%  Kulikowski & Bishop (1981)
%  De Valois, Albrecht, & Thorell (1982)
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2);
maxGridR_deg = max(abs(gridX(:))); %#ok

num_ev_ori  = evS.num_orient; %#ok
num_ev_spf  = evS.num_spFreq;
%num_ev_type = evS.num_cellType; % [1 Complex cell, 4 Simple cells]


%% Set up simulation
%- Get indices of neurons for analysing the results
neuron_orient_deg = specs_simulation.neuron_orient_deg;
neuron_spFreq_cpd = specs_simulation.neuron_spFreq_cpd;
idx_neuron_orient = FindClosestIdx(evS.domain_orient_deg,neuron_orient_deg);
idx_neuron_spFreq = FindClosestIdx(evS.domain_spFreq_cpd,neuron_spFreq_cpd);
idx_neuron_rfLoc = 1;

%- Set specs of stimuli
specs_images.luminance_contrast_Cent = specs_simulation.contrast_pcent/100;
specs_images.radius_deg = specs_simulation.diameter_deg/2;

specs_images.orient_deg = neuron_orient_deg;
specs_images.condition_spFreq_oct = -3:0.05:+3;
specs_images.condition_spFreq_cpd = 2.^(specs_images.condition_spFreq_oct);
numConditions = size(specs_images.condition_spFreq_cpd,2);

% Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast_Cent;
specs_grating.orientation_deg = specs_images.orient_deg;
for c1=1:numConditions
    specs_grating.frequency = specs_images.condition_spFreq_cpd(c1);
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


bandwidth_cpd=NaN(num_models,1);
bandwidth_oct=NaN(num_models,1);

fullH_idx  =NaN(num_models,1);
halfH_idx  =NaN(num_models,2);
fullH_yval =NaN(num_models,1);
halfH_yval =NaN(num_models,1);
fullH_xval_cpd =NaN(num_models,1);
halfH_xval_cpd =NaN(num_models,2);
fullH_xval_oct =NaN(num_models,1);
halfH_xval_oct =NaN(num_models,2);

for tp=1:num_models
    [fullH_yval(tp,1),fullH_idx(tp,1),...
     halfH_yval(tp,1),halfH_idx(tp,:)]=DeriveFWHH(squeeze(results(tp,:))');

    fullH_xval_cpd(tp,1) = ValueFloatIdx(specs_images.condition_spFreq_cpd,fullH_idx(tp,1));
    halfH_xval_cpd(tp,1) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,1));
    halfH_xval_cpd(tp,2) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,2));

    fullH_xval_oct(tp,1) = ValueFloatIdx(specs_images.condition_spFreq_oct,fullH_idx(tp,1));
    halfH_xval_oct(tp,1) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,1));
    halfH_xval_oct(tp,2) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,2));

    bandwidth_cpd(tp) = halfH_xval_cpd(tp,2)-halfH_xval_cpd(tp,1);
    bandwidth_oct(tp) = halfH_xval_oct(tp,2)-halfH_xval_oct(tp,1);
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Spatial frequency';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_spFreq_oct;
retval_mainResults.xlabel = 'Spatial frequency';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'Linear-Rect';'Exponential';'Divisive-Normalization'};

retval_mainResults.bandwidth_cpd = bandwidth_cpd;
retval_mainResults.bandwidth_oct = bandwidth_oct;

retval_mainResults.marking_y = NaN(num_models,3);
retval_mainResults.marking_y(:,1) = fullH_yval;
retval_mainResults.marking_y(:,2) = halfH_yval;
retval_mainResults.marking_y(:,3) = halfH_yval;

retval_mainResults.marking_x(:,1)   = fullH_xval_oct;
retval_mainResults.marking_x(:,2:3) = halfH_xval_oct;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'plot(retval_mainResults.xval,squeeze(retval_mainResults.data(:,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'legend(retval_mainResults.legend_items);'...
                                ...%'hold on;'...
                                ...%'plot(squeeze(retval_mainResults.marking_x(temp_type,:,:)),squeeze(retval_mainResults.marking_y(temp_type,:,:)),''o'');'...
                                ...%'hold off;'...
                                'disp(retval_mainResults.legend_items(:));'...
                                'disp(retval_mainResults.bandwidth_oct(:));'...
                            ];
retval_mainResults.plot = @() eval( retval_mainResults.command );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_SpFreq_CompareModels'' by calling ''retval_mainResults.plot()''.\n'] );

end %%% of file
