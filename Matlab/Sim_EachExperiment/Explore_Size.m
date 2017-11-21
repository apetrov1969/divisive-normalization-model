function [retval_mainResults,...
          retval_allResults] = Explore_Size(M, specs_simulation)

%Explore_Size - Size tuning
%
% Visual stimuli:
%  Luminance gratings with various sizes
%
% Refrerences
%  De Angelis, Freeman, & Ohzawa (1994)
%  Gieselmann & Thiele (2008)
%  Jones, Grieve, Wang, & Sillito (2001)
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

specs_images.orient_deg = neuron_orient_deg;
specs_images.spFreq_cpd = neuron_spFreq_cpd;
specs_images.radius_step_deg = stS.degPerPixel;
specs_images.condition_radius_deg = stS.degPerPixel:specs_images.radius_step_deg:maxGridR_deg;
specs_images.condition_diameter_deg = specs_images.condition_radius_deg*2;
numConditions = size(specs_images.condition_radius_deg,2); % [1,imSize]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast;
specs_grating.orientation_deg = specs_images.orient_deg;
specs_grating.frequency       = specs_images.spFreq_cpd;
for c1=1:numConditions % size
    % Outer radius
    baseGrating = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
    tempImg=baseGrating;
    tempImg(gridR_deg>(specs_images.condition_radius_deg(c1)))=stS.bgLum;
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
results = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,imSize]

%- Find peaks
peak_x=NaN(1,num_ev_type); % [neTyp]
peak_y=NaN(1,num_ev_type); % [neTyp]
for tp=1:num_ev_type
    [peak_y(tp),id] = max(squeeze(results(tp,:,1)));
    if id~=1 && id~=size(specs_images.condition_radius_deg,2)
        peak_x(tp) = specs_images.condition_diameter_deg(id);
    end
end

%- Set 'retval_mainResults'
retval_mainResults.name = 'Size';
retval_mainResults.name = 'Grating size';
retval_mainResults.data = results;
retval_mainResults.xval = specs_images.condition_diameter_deg;
retval_mainResults.xlabel = 'Diameter';
retval_mainResults.ylabel = 'Firing Rate';
retval_mainResults.peak_y = peak_y;
retval_mainResults.peak_x = peak_x;

retval_mainResults.log_name   = 'Log(Size)';
retval_mainResults.log_xlabel = 'Log(Diameter)';

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [  'figure(''Name'',retval_mainResults.name);'...
                                'semilogx(retval_mainResults.xval,squeeze(retval_mainResults.data(temp_type,:)),''-'');'...
                                'title(retval_mainResults.name);'...
                                'xlabel(retval_mainResults.xlabel);'...
                                'ylabel(retval_mainResults.ylabel);'...
                                'hold on;'...
                                'semilogx(retval_mainResults.peak_x(temp_type),'...
                                '     retval_mainResults.peak_y(temp_type), ''o'');'...
                                'hold off;'...
                                '[retval_mainResults.peak_x(temp_type)'...
                                ' retval_mainResults.peak_y(temp_type)]'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Size'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file

