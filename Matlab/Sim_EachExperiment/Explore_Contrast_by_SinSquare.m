function [retval_mainResults,...
          retval_allResults] = Explore_Contrast_by_SinSquare(M,specs_simulation)

%Explore_Contrast_by_SinSquare - Contrasts sensitivity (Sinusoidal/Square)
%
% Visual stimuli:
%  Luminance sinusoidal/square gratings(contrasts) 
%
% References:
%  Pollen & Ronner (1982, 1983): Sinusoidal vs Square
%  Schiller, Finlay, & Volman (1976): Bar vs Dot
%  Morrone, Burr, & Maffei (1982): Bar vs Dot
%  Albrecht, De Valois, & Thorell (1980): Bar vs Grating
%  Watkins & Berkley (1974): Bar vs Grating
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
specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.orient_deg = neuron_orient_deg;
specs_images.condition_wave = 1:2; % Sinusoidal and Square
specs_images.condition_contrast_log10 = 0:0.1:2;
specs_images.condition_contrast = 10.^specs_images.condition_contrast_log10 ./100;
numConditions = [size(specs_images.condition_contrast,2),...
                 size(specs_images.condition_wave,2)]; % [imCon,imSinRec]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imCon,imSinRec]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.frequency       = specs_images.spFreq_cpd;
specs_grating.orientation_deg = specs_images.orient_deg;
for c1=1:numConditions(1) % contrast
    specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast(c1);
    baseGrating = Grating2D(specs_grating,gridX,gridY);
    
    % Sinusoidal
    tempImg = baseGrating+stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg) = stS.bgLum;
    setStimuli(:,:,c1,1) = tempImg;
        
	% Square
    %tempImg = heaviside(baseGrating)-0.5; % 0~+1 -> -0.5~+0.5
    tempImg = baseGrating;
    tempImg(tempImg>0) = 0.5;
    tempImg(tempImg<0) = -0.5;
    tempImg(tempImg==0) = 0;
    
    tempImg = tempImg*2 * specs_grating.amplitude + stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg) = stS.bgLum;
    setStimuli(:,:,c1,2) = tempImg;
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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imCon,imSinRec]

%- Set 'retval_mainResults'
retval_mainResults.name = 'Contrast by Sin/Square';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_contrast;
retval_mainResults.xlabel = 'Contrast';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'sin';'rec'};
retval_mainResults.legend_loc   = 'NorthWest';

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
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Contrast_by_SinSquare'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file

