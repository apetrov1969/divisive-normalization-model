function [retval_mainResults,...
          retval_allResults] = Explore_Simple_Bar(M, specs_simulation)

%Explore_Simple_Bar - Measuring simple-cell receptive field with a bar (position)
%
% Visual stimuli:
%  (1) Single light/dark bars
%  (2) Grating with various spatial-frequencies
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
gridR_deg = sqrt(gridX.^2+gridY.^2);
maxGridR_deg = max(abs(gridX(:)));

num_ev_ori  = evS.num_orient; %#ok
num_ev_spf  = evS.num_spFreq; %#ok
num_ev_type = evS.num_cellType; % [1 Complex cell, 4 Simple cells]

%- Get indices of neurons for analysing the results
neuron_orient_deg = specs_simulation.neuron_orient_deg;
neuron_spFreq_cpd = specs_simulation.neuron_spFreq_cpd;
idx_neuron_orient = FindClosestIdx(evS.domain_orient_deg,neuron_orient_deg);
idx_neuron_spFreq = FindClosestIdx(evS.domain_spFreq_cpd,neuron_spFreq_cpd);
idx_neuron_rfLoc = 1;


%%
%% (1) Bars
%%

%% Set up simulation
%- Set specs of stimuli
specs_images.size_bar_width = specs_simulation.bar_width;

specs_images.back=0.5;
specs_images.light=1;
specs_images.dark=0;
specs_images.x_step = stS.degPerPixel;
specs_images.condition_positionX_deg = -maxGridR_deg:specs_images.x_step:maxGridR_deg;
specs_images.condition_positionY_deg = specs_images.condition_positionX_deg;

numX = size(specs_images.condition_positionX_deg,2);
numConditions = [numX,2];

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
half_width=specs_images.size_bar_width/2;
for px=1:numX % position
    tempImg = NaN([stS.imageSize_pix]);
    tempImg(:)=specs_images.back;
    tempImg(abs(gridX-specs_images.condition_positionX_deg(px)) < half_width) = specs_images.light;
    setStimuli(:,:,px,1) = tempImg;
    
    tempImg = NaN([stS.imageSize_pix]);
    tempImg(:)=specs_images.back;
    tempImg(abs(gridX-specs_images.condition_positionX_deg(px)) < half_width) = specs_images.dark;
    setStimuli(:,:,px,2) = tempImg;
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
    retval_allResults.resp1_CentComp = resp_CentComp;
    retval_allResults.resp1_CentSimp = resp_CentSimp;
    retval_allResults.resp1_SuppChan = resp_SuppChan;
    retval_allResults.resp1_DivNorm  = resp_DivNorm;
    retval_allResults.resp1_SuppDrv  = resp_SuppDrv;
    retval_allResults.resp1_StimDrv  = resp_StimDrv;
end


%% Analyze results
%- Organize results
results1 = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,img]
results1 = reshape(results1,[num_ev_type,numConditions(1),numConditions(2)]); % [neTyp,imPosX,Light/Dark]
results1(:,:,3) = results1(:,:,1)-results1(:,:,2);

%- Set 'retval_mainResults'
retval_mainResults.name1 = 'Line weighting function';
retval_mainResults.data1 = results1;
retval_mainResults.xval1 = specs_images.condition_positionX_deg;
retval_mainResults.xlabel1 = 'Position';
retval_mainResults.ylabel1 = 'Firing Rate';
retval_mainResults.legend_items1 = {'Light';'Dark';'L-D'};
retval_mainResults.legend_loc1   = 'NorthWest';
                          
                          
%%
%% (2) Gratings
%%

clear specs_images specs_grating;
clear setStimuli;
clear numConditions;
clear resp_CentComp resp_CentSimp resp_SuppChan;
clear resp_DivNorm  resp_SuppDrv  resp_StimDrv;

%% Set up simulation
%- Set specs of stimuli
%- Set specs of stimuli
specs_images.luminance_contrast = specs_simulation.contrast_pcent/100;
specs_images.radius_deg = specs_simulation.diameter_deg/2;

specs_images.orient_deg = neuron_orient_deg;
specs_images.condition_spFreq_oct = -3:0.05:+3;
specs_images.condition_spFreq_cpd = 2.^(specs_images.condition_spFreq_oct);
numConditions = size(specs_images.condition_spFreq_cpd,2); % imSpf

%- Generate stimuli
setStimuli = NaN([stS.imageSize_pix, numConditions]);
specs_grating.type = 'Cos';
specs_grating.phase_deg = 0;
specs_grating.amplitude = stS.rangeLum/2 * specs_images.luminance_contrast;
specs_grating.orientation_deg = specs_images.orient_deg;
for c1=1:numConditions
    specs_grating.frequency = specs_images.condition_spFreq_cpd(c1);

    tempImg = Grating2D(specs_grating,gridX,gridY)+stS.bgLum;
    tempImg(gridR_deg>specs_images.radius_deg) = stS.bgLum;

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
    retval_allResults.resp2_CentComp = resp_CentComp;
    retval_allResults.resp2_CentSimp = resp_CentSimp;
    retval_allResults.resp2_SuppChan = resp_SuppChan;
    retval_allResults.resp2_DivNorm  = resp_DivNorm;
    retval_allResults.resp2_SuppDrv  = resp_SuppDrv;
    retval_allResults.resp2_StimDrv  = resp_StimDrv;
end


%% Analyze results
%- Organize results
results2 = squeeze(resp_DivNorm(idx_neuron_orient,idx_neuron_spFreq,idx_neuron_rfLoc,:,:)); % [neOri,neSpf,neRfl,neTyp,img] -> [neTyp,img]

maintained_discharge = evS.divNorm_rateScale_sps * (evS.divNorm_baselineConst^evS.divNorm_exponentNum)/(evS.divNorm_semisaturConst^evS.divNorm_exponentDen);

retval_mainResults.data2a = results2;
retval_mainResults.data2b = results2 - maintained_discharge;
retval_mainResults.xval2 = specs_images.condition_spFreq_oct;
retval_mainResults.maintained_discharge = maintained_discharge;
retval_mainResults.name2   = 'Spatial frequency';
retval_mainResults.xlabel2 = 'Spatial frequency';
retval_mainResults.ylabel2 = 'Firing Rate';


%%
%% (3) Derive frequency tuning curve
%%

measuredMap = squeeze(results1(:,:,3)); % [neTyp,imPosX,Light/Dark/L-D] -> [neTyp,imPosX]
% for tp=1:num_ev_type
%     measuredMap(tp,:) = measuredMap(tp,:)./max(abs(measuredMap(tp,:)));
% end
xval = retval_mainResults.xval1;

gauss1D = @(grid,sig,cen) ( exp(-(grid-cen).^2 ./ (2*sig.^2)) );
gabor_cos = @(sc,sf,sg) (sc .* gauss1D(xval,sg, 0) .* cos(xval*2.*pi.*sf));
gabor_sin = @(sc,sf,sg) (sc .* gauss1D(xval,sg, 0) .* sin(xval*2.*pi.*sf));

sum_sqdiff_2 = @(v) mean((gabor_cos(v(1),v(2),v(3)) - measuredMap(2,:)).^2);
sum_sqdiff_3 = @(v) mean((gabor_sin(v(1),v(2),v(3)) - measuredMap(3,:)).^2);
sum_sqdiff_4 = @(v) mean((gabor_cos(v(1),v(2),v(3)) - measuredMap(4,:)).^2);
sum_sqdiff_5 = @(v) mean((gabor_sin(v(1),v(2),v(3)) - measuredMap(5,:)).^2);

temp_sg = M.EarlyVis_filters.sigma_clrf_width(idx_neuron_spFreq);
temp_sf = neuron_spFreq_cpd;
opt_values=cell(1,5);

tp=2; temp_sc = max(abs(measuredMap(tp,:))); temp_values = [temp_sc,temp_sf,temp_sg]; opt_values{tp} = fminsearch(sum_sqdiff_2,temp_values);
tp=3; temp_sc =-max(abs(measuredMap(tp,:))); temp_values = [temp_sc,temp_sf,temp_sg]; opt_values{tp} = fminsearch(sum_sqdiff_3,temp_values);
tp=4; temp_sc =-max(abs(measuredMap(tp,:))); temp_values = [temp_sc,temp_sf,temp_sg]; opt_values{tp} = fminsearch(sum_sqdiff_4,temp_values);
tp=5; temp_sc = max(abs(measuredMap(tp,:))); temp_values = [temp_sc,temp_sf,temp_sg]; opt_values{tp} = fminsearch(sum_sqdiff_5,temp_values);

scaleFctr=NaN(1,5);
cent_cpd =NaN(1,5);
sigma_deg=NaN(1,5);
sigma_cpd=NaN(1,5);
for tp=2:num_ev_type
    scaleFctr(tp) =opt_values{tp}(1);
    cent_cpd(tp)  =opt_values{tp}(2);
    sigma_deg(tp) =opt_values{tp}(3);
    sigma_cpd(tp) =(1/(2*pi)) ./ sigma_deg(tp);
end

derived_curves=NaN(num_ev_type,size(specs_images.condition_spFreq_cpd,2));
for tp=2:num_ev_type
    derived_curves(tp,:) = gauss1D(specs_images.condition_spFreq_cpd,sigma_cpd(tp),cent_cpd(tp));
end
measured_curves=NaN(num_ev_type,size(specs_images.condition_spFreq_cpd,2),2);
measured_curves(:,:,1) = retval_mainResults.data2a;
measured_curves(:,:,2) = retval_mainResults.data2b;

%- Compute bandwidths of derived/measured spf tuning curves
num_curves = 3;

bandwidth_cpd=NaN(num_ev_type,num_curves); % measured/measured-mdc/derived
bandwidth_oct=NaN(num_ev_type,num_curves);

fullH_idx  =NaN(num_ev_type,1,num_curves);
halfH_idx  =NaN(num_ev_type,2,num_curves);
fullH_yval =NaN(num_ev_type,1,num_curves);
halfH_yval =NaN(num_ev_type,1,num_curves);
fullH_xval_cpd =NaN(num_ev_type,1,num_curves);
halfH_xval_cpd =NaN(num_ev_type,2,num_curves);
fullH_xval_oct =NaN(num_ev_type,1,num_curves);
halfH_xval_oct =NaN(num_ev_type,2,num_curves);

for tp=1:num_ev_type
    cv=1;
    [fullH_yval(tp,1,cv),fullH_idx(tp,1,cv),...
     halfH_yval(tp,1,cv),halfH_idx(tp,:,cv)]=DeriveFWHH(squeeze(measured_curves(tp,:,cv)));
    cv=2;
    [fullH_yval(tp,1,cv),fullH_idx(tp,1,cv),...
     halfH_yval(tp,1,cv),halfH_idx(tp,:,cv)]=DeriveFWHH(squeeze(measured_curves(tp,:,cv)));
    cv=3;
    [fullH_yval(tp,1,cv),fullH_idx(tp,1,cv),...
     halfH_yval(tp,1,cv),halfH_idx(tp,:,cv)]=DeriveFWHH(squeeze(derived_curves(tp,:)));
 
    for cv=1:num_curves
        fullH_xval_cpd(tp,1,cv) = ValueFloatIdx(specs_images.condition_spFreq_cpd,fullH_idx(tp,1,cv));
        halfH_xval_cpd(tp,1,cv) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,1,cv));
        halfH_xval_cpd(tp,2,cv) = ValueFloatIdx(specs_images.condition_spFreq_cpd,halfH_idx(tp,2,cv));

        fullH_xval_oct(tp,1,cv) = ValueFloatIdx(specs_images.condition_spFreq_oct,fullH_idx(tp,1,cv));
        halfH_xval_oct(tp,1,cv) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,1,cv));
        halfH_xval_oct(tp,2,cv) = ValueFloatIdx(specs_images.condition_spFreq_oct,halfH_idx(tp,2,cv));

        bandwidth_cpd(tp,cv) = halfH_xval_cpd(tp,2,cv)-halfH_xval_cpd(tp,1,cv);
        bandwidth_oct(tp,cv) = halfH_xval_oct(tp,2,cv)-halfH_xval_oct(tp,1,cv);
    end
end

retval_mainResults.curve1 = derived_curves;
retval_mainResults.legend_itemsA = {'Measured';'Measured-mdc';'Derived'};
retval_mainResults.legend_locA   = 'NorthWest';
retval_mainResults.scaleFctr = scaleFctr;
retval_mainResults.cent_cpd  = cent_cpd;
retval_mainResults.sigma_deg = sigma_deg;
retval_mainResults.sigma_cpd = sigma_cpd;

retval_mainResults.bandwidth_cpd = bandwidth_cpd;
retval_mainResults.bandwidth_oct = bandwidth_oct;

retval_mainResults.markingA_y = NaN(num_ev_type,3,num_curves);
retval_mainResults.markingA_y(:,1,:) = fullH_yval;
retval_mainResults.markingA_y(:,2,:) = halfH_yval;
retval_mainResults.markingA_y(:,3,:) = halfH_yval;

retval_mainResults.markingA_x(:,1,:)   = fullH_xval_oct;
retval_mainResults.markingA_x(:,2:3,:) = halfH_xval_oct;


%% Plot A
complex = 1; %#ok
simple = @(ph) (mod(round(ph/90),4)+2); %#ok

retval_mainResults.commandA = [
                                'figure(''Name'',retval_mainResults.name2);'...
                                'title(retval_mainResults.name2);'...
                                'xlabel(retval_mainResults.xlabel2);'...
                                'ylabel(retval_mainResults.ylabel2);'...
                                'hold on;'...
                                'plot(retval_mainResults.xval2,squeeze(retval_mainResults.data2a(temp_type,:)./max(retval_mainResults.data2a(temp_type,:))),''-b'');'...
                                'plot(retval_mainResults.xval2,squeeze(retval_mainResults.data2b(temp_type,:)./max(retval_mainResults.data2b(temp_type,:))),''-g'');'...
                                'plot(retval_mainResults.xval2,squeeze(retval_mainResults.curve1(temp_type,:)./max(retval_mainResults.curve1(temp_type,:))),''-r'');'...
                                'legend(retval_mainResults.legend_itemsA);'...
                                'legend(''Location'',retval_mainResults.legend_locA);'...
                                'plot(squeeze(retval_mainResults.markingA_x(temp_type,:,1)),  squeeze(retval_mainResults.markingA_y(temp_type,:,1))./max(retval_mainResults.data2a(temp_type,:)),''ob'');'...
                                'plot(squeeze(retval_mainResults.markingA_x(temp_type,:,2)),  squeeze(retval_mainResults.markingA_y(temp_type,:,2))./max(retval_mainResults.data2b(temp_type,:)),''og'');'...
                                'plot(squeeze(retval_mainResults.markingA_x(temp_type,:,3)),  squeeze(retval_mainResults.markingA_y(temp_type,:,3))./max(retval_mainResults.curve1(temp_type,:)),''or'');'...
                                'hold off;'...
                                'fprintf(''Bandwidths(Measured, Measured-mdc, Derived)=\n'');'...
                                'disp(retval_mainResults.bandwidth_oct(temp_type,:));'...
                             ];
retval_mainResults.plotA = @(tp) eval( strrep( retval_mainResults.commandA, 'temp_type', num2str(eval(strcat(tp))) ) );


%% Plot B
retval_mainResults.legend_itemsB = {'Measured';'Light-Dark';'Fitted'};
retval_mainResults.legend_locB   = 'NorthWest';
                            
retval_mainResults.commandB = [ 
                                'figure(''Name'',retval_mainResults.name1);'...
                                'title(retval_mainResults.name1);'...
                                'xlabel(retval_mainResults.xlabel1);'...
                                'ylabel(retval_mainResults.ylabel1);'...
                                'hold on;'...
                                'plot(retval_mainResults.xval1,squeeze(retval_mainResults.data1(temp_type,:,1)),''-b'');'...
                                'plot(retval_mainResults.xval1,squeeze(retval_mainResults.data1(temp_type,:,3)),''-g'');'...
                                'if(temp_type==2 || temp_type==4) '...
                                    'plot(retval_mainResults.xval1,gabor_cos(retval_mainResults.scaleFctr(temp_type),'...
                                                                            'retval_mainResults.cent_cpd (temp_type),'...
                                                                            'retval_mainResults.sigma_deg(temp_type)),''-r'');'...
                                'else '...
                                    'plot(retval_mainResults.xval1,gabor_sin(retval_mainResults.scaleFctr(temp_type),'...
                                                                            'retval_mainResults.cent_cpd (temp_type),'...
                                                                            'retval_mainResults.sigma_deg(temp_type)),''-r'');'...
                                'end;'...
                                'legend(retval_mainResults.legend_itemsB);'...
                                'legend(''Location'',retval_mainResults.legend_locB);'...
                                'plot(retval_mainResults.xval1,-squeeze(retval_mainResults.data1(temp_type,:,2)),''-b'');'...
                                'hold off;'...
                             ];
retval_mainResults.plotB = @(tp) eval( strrep( retval_mainResults.commandB, 'temp_type', num2str(eval(strcat(tp))) ) );

%- Set 'comment'
retval_mainResults.comment = sprintf( ['Plot the results of ''Explore_Simple_Bar'' by calling ''retval_mainResults.plotA/B(tp)''.\n'...
                                       'where ''tp'' represent a type of neurons: ''complex''/''simple(0,90,180 or 270)''\n'] );

end %%% of file