function [retval_resp_divNorm,...
          retval_resp_suppDrv,...
          retval_resp_stimDrv] = DMPL_EarlyVis_MeanChannelResp(M,resp_centComp,...
                                                                 resp_centSimp,...
                                                                 resp_suppChan)


%%  Retrieve some oft-used parameters
num_chan_SpFreq = M.EarlyVis_spec.num_channel_spFreq; % chanSpf
num_orient = M.EarlyVis_spec.num_orient; % neurOri
num_spFreq = M.EarlyVis_spec.num_spFreq; % neurSpf
num_rfLoc  = M.EarlyVis_spec.num_rfLoc;  % neurRfl
num_types  = M.EarlyVis_spec.num_cellType; % 1 Complex cell + 4 Simple cells(phase:0,90,180,270)
num_images = size(resp_suppChan,5);


%% Calibrattion 1/2: Individual channels
resp_centSimp = resp_centSimp ./ M.EarlyVis_calibration.calibrFactor_numerator_simp; % [neurOri, neurSpf, neurRfl, neurPhs, img]
resp_centComp = resp_centComp ./ M.EarlyVis_calibration.calibrFactor_numerator_comp; % [neurOri, neurSpf, neurRfl, img]


%% Retrieve parameters
divNorm_exponentNum = M.EarlyVis_spec.divNorm_exponentNum;
divNorm_exponentDen = M.EarlyVis_spec.divNorm_exponentDen;
divNorm_baselineConst  = M.EarlyVis_spec.divNorm_baselineConst;
divNorm_semisaturConst = M.EarlyVis_spec.divNorm_semisaturConst;
divNorm_rateScale_sps  = M.EarlyVis_spec.divNorm_rateScale_sps;
orient_kernel_weightMatrix = M.EarlyVis_divNorm.orient_kernel_weightMatrix;
spFreq_kernel_weightMatrix = M.EarlyVis_divNorm.spFreq_kernel_weightMatrix;


%% Compute Stimulus drive of Divisive-normalization
retval_resp_stimDrv = NaN(num_spFreq,num_orient,num_rfLoc,num_types,num_images); % [neurSpf,neurOri,neurRfl,neurCompSimp,img]

%- Complex cell: channels_Comp[neurOri, neurSpf, neurRfl, img]
temp_chanComp = permute(resp_centComp,[2,1,3,4]);       % [neurSpf,neurOri,neurRfl,img]
temp_chanComp = reshape(temp_chanComp,[num_spFreq,num_orient,num_rfLoc,1,num_images]); % [neurSpf,neurOri,neurRfl,1,img]
retval_resp_stimDrv(:,:,:,1,:) = temp_chanComp(:,:,:,:,:);   % [neurSpf,neurOri,neurRfl,neurCompSimp,img]

%- Simple cell: channels_Simp[neurOri, neurSpf, neurRfl, neurPhs, img]
temp_chanSimp = permute(resp_centSimp,[2,1,3,4,5]);                 % [neurSpf,neurOri,neurRfl,nePhs,img]
retval_resp_stimDrv(:,:,:,2:num_types,:) = temp_chanSimp(:,:,:,:,:);   % [neurSpf,neurOri,neurRfl,neurCompSimp,img]

               
%% Compute Suppressive drive of Divisive-normalization    
temp_suppDrive = resp_suppChan; % [chanOri, chanSpf, neurSpf, neurRfl, img]

%- Weighted sum across the orientation channels
temp_suppDrive = reshape(temp_suppDrive,[num_orient,num_chan_SpFreq*num_spFreq*num_rfLoc*num_images]); % [chanOri,chanSpf*neurSpf*neurRfl*img]
temp_suppDrive = orient_kernel_weightMatrix * temp_suppDrive; % [neurOri,chanOri]*[chanOri,chanSpf*neurSpf*neurRfl*img] -> [neurOri,chanSpf*neurSpf*neurRfl*img]

%- Weighted sum across the spatial-frequency channels
temp_suppDrive = reshape(temp_suppDrive,[num_orient,num_chan_SpFreq,num_spFreq,num_rfLoc,num_images]); % [neurOri,chanSpf,neurSpf,neurRfl,img]
temp_suppDrive = permute(temp_suppDrive,[2,1,3,4,5]);                                                  % [chanSpf,neurOri,neurSpf,neurRfl,img]
temp_suppDrive = reshape(temp_suppDrive,[num_chan_SpFreq,num_orient*num_spFreq*num_rfLoc*num_images]); % [chanSpf,neurOri*neurSpf*neurRfl*img]
temp_suppDrive = spFreq_kernel_weightMatrix * temp_suppDrive; % [neurSpf,chanSpf]*[chanSpf,neurOri*neurSpf*neurRfl*img] -> [neurSpf,neurOri*neurSpf*neurRfl*img]

temp_suppDrive = reshape(temp_suppDrive,[num_spFreq,num_orient,num_spFreq,num_rfLoc,num_images]); % [neurSpf,neurOri,neurSpf,neurRfl,img]
temp_suppDrive = permute(temp_suppDrive, [1,3,2,4,5]);                                            % [neurSpf,neurSpf,chanOri,neurRfl,img]

retval_resp_suppDrv = NaN(num_spFreq,num_orient,num_rfLoc,num_images);                               % [neurSpf,neurOri,neurRfl,img]
for sf=1:num_spFreq
    retval_resp_suppDrv(sf,:,:,:) = temp_suppDrive(sf,sf,:,:,:);                                     % [neurSpf,neurOri,neurRfl,img]
end


%% Calibrattion 2/2: Denominator
retval_resp_suppDrv = retval_resp_suppDrv ./ M.EarlyVis_calibration.calibrFactor_denominator;


%% Compute Numerator/Denominator of Divisive-normalization
temp_Numerator = retval_resp_stimDrv + divNorm_baselineConst;	% [neurSpf,neurOri,neurRfl,neurCompSimp,img]
temp_Numerator = (abs(temp_Numerator)+temp_Numerator)./2;	% [neurSpf,neurOri,neurRfl,neurCompSimp,img]
divNormNumerator   = temp_Numerator.^divNorm_exponentNum;     % [neurSpf,neurOri,neurRfl,neurCompSimp,img]
divNormDenominator = retval_resp_suppDrv + divNorm_semisaturConst^divNorm_exponentDen; % [neurSpf,neurOri,neurRfl,img]


%% Calculate and return the normalized ratio
retval_resp_divNorm = NaN(num_spFreq,num_orient,num_rfLoc,num_types,num_images);             % [neurSpf,neurOri,neurRfl,neurCompSimp,img]
tempDenominator = reshape(divNormDenominator,[num_spFreq,num_orient,num_rfLoc,1,num_images]);   % [neurSpf,neurOri,neurRfl,1,img]
for tp=1:num_types
    retval_resp_divNorm(:,:,:,tp,:) = divNormNumerator(:,:,:,tp,:)./tempDenominator;
end
retval_resp_divNorm = retval_resp_divNorm.*divNorm_rateScale_sps; % [neurSpf,neurOri,neurRfl,neurCompSimp,img]
retval_resp_divNorm = permute(retval_resp_divNorm,[2,1,3,4,5]);   % [neurOri,neurSpf,neurRfl,neurCompSimp,img]

retval_resp_suppDrv = permute(retval_resp_suppDrv,[2,1,3,4]);     % [neurOri,neurSpf,neurRfl,img]
retval_resp_stimDrv = permute(retval_resp_stimDrv,[2,1,3,4,5]);   % [neurOri,neurSpf,neurRfl,neurCompSimp,img]

%%% Return retval_meanChannelResp and optionally suppressDrive
end %%% of file
