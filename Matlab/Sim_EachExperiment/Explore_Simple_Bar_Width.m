function [retval_mainResults,...
          retval_allResults] = Explore_Simple_Bar_Width(M, specs_simulation)

%Explore_Simple_Bar_Width - Measuring simple-cell receptive field with a bar (position x width)
%
% Visual stimuli:
%  single light/dark bars
%
% Refrerences:
%  TadmorTolhurst (1998)
%  AndrewsPollen (1979)
%  Kulikowski & Vidyasagar (1986)
%  Kulikowski, Marcelja, & Bishop (1982)
%  Kulikowski & Bishop (1981)
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
specs_images.back=0.5;
specs_images.light=1;
specs_images.dark=0;

x_step = stS.degPerPixel;
maxW = 0.5/neuron_spFreq_cpd;
maxPosX = maxGridR_deg - maxW/2;
specs_images.x_step = x_step;
specs_images.condition_barWidth = x_step:x_step:maxW;
specs_images.condition_positionX_deg = -maxPosX:x_step:+maxPosX;

numConditions = [size(specs_images.condition_positionX_deg,2),...
                 size(specs_images.condition_barWidth,2), 2]; % [imPosX,imBarW,Light/Dark]

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
for px=1:numConditions(1) % position
    for wd=1:numConditions(2) % width
        half_width=specs_images.condition_barWidth(wd)/2;

        tempImg = NaN([stS.imageSize_pix]);
        tempImg(:)=specs_images.back;
        tempImg(abs(gridX-specs_images.condition_positionX_deg(px)) < half_width) = specs_images.light;
        setStimuli(:,:,px,wd,1) = tempImg;

        tempImg = NaN([stS.imageSize_pix]);
        tempImg(:)=specs_images.back;
        tempImg(abs(gridX-specs_images.condition_positionX_deg(px)) < half_width) = specs_images.dark;
        setStimuli(:,:,px,wd,2) = tempImg;
    end
end
setStimuli = reshape(setStimuli,[stS.imageSize_pix, numConditions(1)*numConditions(2)*numConditions(3)]); % [y,x,img]


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
results = reshape(results,[num_ev_type,numConditions(1),numConditions(2),numConditions(3)]); % [neTyp,imPosX,imBarW,Light/Dark]
results(:,:,:,3) = results(:,:,:,1)-results(:,:,:,2); % [neTyp,imPosX,imBarW,Light/Dark/L-D]


%% Derive frequency tuning curve

derived_fwhh_oct = NaN(num_ev_type,numConditions(2)); % [neTyp,imBarW]
derived_cent_cpd = NaN(num_ev_type,numConditions(2)); % [neTyp,imBarW]
scaleFctr = NaN(1,num_ev_type);
cent_cpd  = NaN(1,num_ev_type);
sigma_deg = NaN(1,num_ev_type);
sigma_cpd = NaN(1,num_ev_type);
sum_sqdiff = cell(1,5);
for wd=1:numConditions(2) % width
    temp_results = squeeze(results(:,:,wd,3)); % [neTyp,imPosX,imBarW,Light/Dark/L-D] -> [neTyp,imPosX]
    xval = specs_images.condition_positionX_deg;

    gauss1D = @(grid,sig,cen) ( exp(-(grid-cen).^2 ./ (2*sig.^2)) );
    gabor_cos = @(sc,sf,sg) (sc .* gauss1D(xval,sg, 0) .* cos(xval*2.*pi.*sf)); % (scale, spf, sigma)
    gabor_sin = @(sc,sf,sg) (sc .* gauss1D(xval,sg, 0) .* sin(xval*2.*pi.*sf)); % (scale, spf, sigma)

    sum_sqdiff{2} = @(v) mean((gabor_cos(v(1),v(2),v(3)) - temp_results(2,:)).^2);
    sum_sqdiff{3} = @(v) mean((gabor_sin(v(1),v(2),v(3)) - temp_results(3,:)).^2);
    sum_sqdiff{4} = @(v) mean((gabor_cos(v(1),v(2),v(3)) - temp_results(4,:)).^2);
    sum_sqdiff{5} = @(v) mean((gabor_sin(v(1),v(2),v(3)) - temp_results(5,:)).^2);

    temp_sg = M.EarlyVis_filters.sigma_clrf_width(idx_neuron_spFreq);
    temp_sf = neuron_spFreq_cpd;
    opt_values=cell(1,num_ev_type);

    scaleFctr(:) = 0;
    cent_cpd (:) = 0;
    sigma_deg(:) = 0;
    sigma_cpd(:) = 0;
    for tp=2:num_ev_type
    	temp_sc = max(abs(temp_results(tp,:)));
        temp_values = [temp_sc,temp_sf,temp_sg];
        opt_values{tp} = fminsearch(sum_sqdiff{tp},temp_values);
        
        scaleFctr(tp) =opt_values{tp}(1);
        cent_cpd(tp)  =opt_values{tp}(2);
        sigma_deg(tp) =opt_values{tp}(3);
        sigma_cpd(tp) =(1/(2*pi)) ./ sigma_deg(tp);
    end
    
    fwhh_cpd = NormalDistrib_StdDev_to_FWHH(sigma_cpd);
    hh_left  = cent_cpd - fwhh_cpd./2;
    hh_right = cent_cpd + fwhh_cpd./2;
    hh_left  = log2(hh_left);
    hh_right = log2(hh_right);
    derived_fwhh_oct(:,wd) = hh_right-hh_left;
    derived_cent_cpd(:,wd) = cent_cpd;
end


%% Set 'retval_mainResults'
retval_mainResults.name1 = 'Bar position x width';
retval_mainResults.data = results;
retval_mainResults.xval1 = specs_images.condition_positionX_deg;
retval_mainResults.xlabel1 = 'Position';
retval_mainResults.ylabel1 = 'Firing Rate';
retval_mainResults.legend_items = cellstr(num2str(specs_images.condition_barWidth','%.2f°'));
retval_mainResults.legend_loc   = 'NorthWest';

retval_mainResults.name2 = 'Bandwidth (FWHH) x Image probe';
retval_mainResults.xval2 = specs_images.condition_barWidth;
retval_mainResults.xlabel2 = 'Bar width (deg)';
retval_mainResults.ylabel2 = 'Derived bandwidth (oct)';
retval_mainResults.fwhh_oct = derived_fwhh_oct;

retval_mainResults.name3 = 'Preferred frequency x Image probe';
retval_mainResults.xval3 = specs_images.condition_barWidth;
retval_mainResults.xlabel3 = 'Bar width (deg)';
retval_mainResults.ylabel3 = 'Derived frequency (cpd)';
retval_mainResults.cent_cpd = derived_cent_cpd;

%- Set 'plot' function
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok
retval_mainResults.command = [
                                'figure(''Name'',retval_mainResults.name1);'...
                                'plot(retval_mainResults.xval1,squeeze(retval_mainResults.data(temp_type,:,:,3)),''-'');'...
                                'title(retval_mainResults.name1);'...
                                'xlabel(retval_mainResults.xlabel1);'...
                                'ylabel(retval_mainResults.ylabel1);'...
                                'legend(retval_mainResults.legend_items);'...
                                'legend(''Location'',retval_mainResults.legend_loc);'...
                                ...
                                'figure(''Name'',retval_mainResults.name2);'...
                                'plot(retval_mainResults.xval2,squeeze(retval_mainResults.fwhh_oct(temp_type,:)),''-'');'...
                                'title(retval_mainResults.name2);'...
                                'xlabel(retval_mainResults.xlabel2);'...
                                'ylabel(retval_mainResults.ylabel2);'...
                                ...
                                'figure(''Name'',retval_mainResults.name3);'...
                                'plot(retval_mainResults.xval3,squeeze(retval_mainResults.cent_cpd(temp_type,:)),''-'');'...
                                'title(retval_mainResults.name3);'...
                                'xlabel(retval_mainResults.xlabel3);'...
                                'ylabel(retval_mainResults.ylabel3);'...
                            ];
retval_mainResults.plot = @(tp) eval( strrep( retval_mainResults.command, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Simple_Bar_Width'' by calling ''retval_mainResults.plot(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file
