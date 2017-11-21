function [retval_mainResults,...
          retval_allResults] = Explore_PlaidContrast(M, specs_simulation)

%Explore_PlaidContrast - Contrast sensitivity (Plaid)
%
% Visual stimuli:
%  Plaid (Two gratings superimposed) with various contrasts
%
% References:
%  Koch, Jin, Alonso, & Zaidi (2016)
%  Priebe & Ferster (2006)
%  Busse, Wade, & Carandini (2009)
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

specs_images.condition_contrast = 10.^(0:0.1:log10(50))./100;
specs_images.condition_contrast = [specs_images.condition_contrast, 0.5];
specs_images.condition_contrast = sort(specs_images.condition_contrast);

numConditions = size(specs_images.condition_contrast,2); % imCon

% Generate stimuli
setStimuli = NaN([stS.imageSize_pix, 2, numConditions]); % [x,y,imCon]
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.frequency = specs_images.spFreq_cpd;

for c1=1:numConditions % contrast
    specs_grating.amplitude = stS.rangeLum/2 * specs_images.condition_contrast(c1);
    
    specs_grating.orientation_deg = specs_images.orient1_deg;
    grating1 = Grating2D(specs_grating,gridX,gridY);
    specs_grating.orientation_deg = specs_images.orient2_deg;
    grating2 = Grating2D(specs_grating,gridX,gridY);
        
    tempImg = grating1+stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
    setStimuli(:,:,1,c1) = tempImg;
    
    tempImg = grating1+grating2+stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg)=stS.bgLum;
    setStimuli(:,:,2,c1) = tempImg;
    
end % for c1=1:numConditions % contrast
setStimuli = reshape(setStimuli,[stS.imageSize_pix, 2*numConditions]); % [y,x,img]


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
%resp_DivNorm = resp_StimDrv;
results = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,img]

%results = squeeze(resp_SuppDrv(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:)); % [neOri,neSpf,neRfl,img] -> [1,img]_
%results = repmat(results',[5,1]);

results = reshape(results,[num_ev_type,2,numConditions]); % [neTyp,2,imCon]

%- Set 'retval_mainResults'
retval_mainResults.name = 'Plaid contrast';
retval_mainResults.data = results; % [neTyp,imGratingPlaid,imCon]
retval_mainResults.supp_index = 1 - squeeze(results(:,2,:)./results(:,1,:)); % [neTyp,imCon]
retval_mainResults.xval = specs_images.condition_contrast;
retval_mainResults.xlabel = 'Plaid contrast';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.legend_items = {'Grating';'Plaid'};
retval_mainResults.legend_loc   = 'NorthWest';

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  ...'figure(''Name'',retval_mainResults.name);'...
                                ...'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:,:)),''-'');'...
                                ...'xlabel(retval_mainResults.xlabel);'...
                                ...'ylabel(retval_mainResults.ylabel);'...
                                ...'axis([0.01 1 0 inf]);'...
                                ...'legend(retval_mainResults.legend_items);'...
                                ...'legend(''Location'',retval_mainResults.legend_loc);'...
                                ...
                                'figure(''Name'',''Suppression Index'');'...
                                'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.supp_index(temp_type,:)),''-'');'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(''Suppression Index'');'...
                                'axis([0.01 1 -0.1 1]);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_PlaidContrast'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file