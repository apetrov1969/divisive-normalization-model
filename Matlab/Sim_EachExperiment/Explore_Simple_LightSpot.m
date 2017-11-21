function [retval_mainResults,...
          retval_allResults] = Explore_Simple_LightSpot(M, specs_simulation)

%Explore_Simple_LightSpot - Measuring simple-cell receptive field using a light spot
%
% Visual stimuli:
%  Single light spots (see Hubel & Wiesel, 1959)
%
% Refrerences
%  Hubel& Wiesel (1959)
%

%% Retrieve parameters
stS = M.stim_spec;
evS = M.EarlyVis_spec;
gridX = stS.gridX_deg;
gridY = stS.gridY_deg;
gridR_deg = sqrt(gridX.^2+gridY.^2); %#ok
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
specs_images.size_spot_deg = stS.degPerPixel;
specs_images.light=1;
specs_images.dark=0;
specs_images.back=0.5;
specs_images.yx_step = stS.degPerPixel*1;
specs_images.condition_positionX_deg = -maxGridR_deg:specs_images.yx_step:maxGridR_deg;
specs_images.condition_positionY_deg = specs_images.condition_positionX_deg;

numY = size(specs_images.condition_positionY_deg,2);
numX = size(specs_images.condition_positionX_deg,2);
numConditions = [numY,numX];

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
square_half_spot=(specs_images.size_spot_deg/2).^2;
for py=1:numY % position
    for px=1:numX % position
        tempImg = NaN([stS.imageSize_pix]);
        tempImg(:)=specs_images.back;
        tempImg((gridY-specs_images.condition_positionY_deg(py)).^2 ...
               +(gridX-specs_images.condition_positionX_deg(px)).^2 < square_half_spot) = specs_images.light;
        setStimuli(:,:,py,px) = tempImg;
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
results = permute(results,[2,1]); % [img,neTyp]
results = reshape(results,[numY,numX,num_ev_type]); % [posY,posX,neTyp]

%- Set color map
cm1=[0:255; 0:255; 255*ones(1,256)];
cm2=[255*ones(1,256); 255:-1:0; 255:-1:0];
cm3=[cm1';cm2'];
cm3=cm3/255;

%- Set 'retval_mainResults'
retval_mainResults.name = 'Receptive field (single spots)';
retval_mainResults.data = results;
retval_mainResults.colormap = cm3;
retval_mainResults.maxResp = max(results(:));
retval_mainResults.minResp = min(results(:));
retval_mainResults.bgResp = squeeze(results(1,1,:));
retval_mainResults.coordinatesStyle = stS.coordinatesStyle; % 'ij' or 'xy'

%- Set 'plot' function
tp=0; %#ok
retval_mainResults.command = [  'for tp=1:num_ev_type;'...
                                    'figure(''Name'',retval_mainResults.name);'...
                                    'imagesc(squeeze(retval_mainResults.data(:,:,tp)));'...
                                    'set(gca,''xtick'',[],''ytick'',[]);'...
                                    'axis(retval_mainResults.coordinatesStyle);'...
                                    'axis(''image'');'...
                                    'colorbar;'...
                                    'colormap(retval_mainResults.colormap);'...
                                    'truesize;'...
                                'end;'...
                            ];
retval_mainResults.plot = @(tp) eval( retval_mainResults.command );

%- Set 'comment'
retval_mainResults.comment = sprintf( 'Plot the results of ''Explore_Simple_LightSpot'' by calling ''retval_mainResults.plot()''.' );

end %%% of file


