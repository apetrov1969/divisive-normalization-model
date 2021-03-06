function [retval_mainResults,...
          retval_allResults] = Explore_SurrSpFreq(M, specs_simulation)

%Explore_SurrSpFreq - Surr spatial-frequency tuning
%
% Visual stimuli:
%  Luminance gratings with various frequencies in center/surround
%
% References:
%  Maffei & Fiorentini (1976)
%  Li & Li (1994)
%  De Angelis, Freeman, & Ohzawa (1994)
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
specs_images.luminance_contrast = specs_simulation.contrast_pcent/100;
specs_images.radius_Cent_deg = specs_simulation.center_diameter_deg/2;
specs_images.radius_inner_Surr_deg = specs_simulation.surround_inner_diameter_out_deg/2;
specs_images.radius_outer_Surr_deg = specs_simulation.surround_outer_diameter_out_deg/2;

specs_images.spFreq_base_cpd = neuron_spFreq_cpd;
specs_images.orient_deg = neuron_orient_deg;

specs_images.condition_spFreq_oct = -3:0.05:+3;
specs_images.condition_spFreq_cpd = 2.^(specs_images.condition_spFreq_oct);
numConditions = [size(specs_images.condition_spFreq_cpd,2),2]; % [imSpf,cent/surr]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions(1), numConditions(2)]);
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast;
specs_grating.frequency = specs_images.spFreq_base_cpd;
specs_grating.orientation_deg = specs_images.orient_deg;
baseGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
for c=1:numConditions(1) % spatial frequency
    specs_grating.frequency = specs_images.condition_spFreq_cpd(c);
    specs_grating.phase_deg = 180;
    varGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
    
    % Center(varGrating)
    tempImg=varGrating;
    tempImg(gridR_deg>(specs_images.radius_Cent_deg))=stS.bgLum;
    setStimuli(:,:,c,1) = tempImg;
        
    % Center(baseGrating) + Surround(varGrating)
    tempImg=varGrating;
    tempImg(gridR_deg>specs_images.radius_outer_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_inner_Surr_deg)=stS.bgLum;
    tempImg(gridR_deg<specs_images.radius_Cent_deg)=baseGrating(gridR_deg<(specs_images.radius_Cent_deg));
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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imSpf,cent/surr]

%- Set 'retval_mainResults'
retval_mainResults.name = 'Surround suppression (Frequency)';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_spFreq_oct;
retval_mainResults.xlabel = 'Spatial frequency';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'Cent only';'Cent+Surr'};
retval_mainResults.legend_loc   = 'NorthWest';

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
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_SurrSpFreq'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
