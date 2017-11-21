function [retval_mainResults,...
          retval_allResults] = Explore_Contrast(M,specs_simulation)

%Explore_Contrast - Contrasts sensitivity
%
% Visual stimuli:
%  Luminance gratings(contrasts)
%
% References:
%  Carandini (2004)
%  Carandini, Heeger, & Movshon (1997)
%  Albrecht, Geisler, & Crane (2003)
%  Albrecht, Hamilton (1982)
%  Itti, Koch, & Braun (2000)
%  Tolhurst & Dean (1987)
%  Schumer & Movshon (1984)
%  Li & Creutzfeldt (1984)
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
num_ev_type = evS.num_cellType;  %#ok % [1 Complex cell, 4 Simple cells]

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
specs_images.condition_contrast_log10 = 0:0.1:2;
specs_images.condition_contrast   = 10.^specs_images.condition_contrast_log10 ./100;
numConditions = size(specs_images.condition_contrast,2); % [1,imCon]

% Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]); % [x,y,imCon]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.orientation_deg = specs_images.orient_deg;
specs_grating.frequency       = specs_images.spFreq_cpd;
for c1=1:numConditions % contrast
    specs_grating.amplitude = stS.rangeLum/2 .* specs_images.condition_contrast(c1);
    tempImg = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
    tempImg(gridR_deg > specs_images.radius_deg) = stS.bgLum;
    setStimuli(:,:,c1) = tempImg;
end
setStimuli = reshape(setStimuli,[stS.imageSize_pix, numConditions]); % [y,x,img]


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


%% Analyze results : retval_mainResults
%- Organize results
results = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,imCon]

%- Contrast x Orientation
retval_mainResults.name = 'Contrast';
retval_mainResults.data = results;
retval_mainResults.xval1 = specs_images.condition_contrast;
retval_mainResults.xlabel1 = 'Log10Contrast';
retval_mainResults.xval2 = specs_images.condition_contrast.*100;
retval_mainResults.xlabel2 = 'Contrast (%)';
retval_mainResults.ylabel = 'Firing Rate';

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'semilogx(retval_mainResults.xval1,squeeze(retval_mainResults.data(temp_type,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel1);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'axis([0 2 0 50]);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Contrast'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
